"""Orchestration helpers for the GenBank to BioDataset pipeline."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from .cache import get_cached_text, set_cached_text
from .featurizer import clean_sequence, compute_features
from .ncbi_client import NCBIClient, NCBIClientConfig
from .exporter import write_csv, write_fasta, write_jsonl


@dataclass(slots=True)
class RunConfig:
    query: str
    limit: int
    email: str
    tool: str
    api_key: str | None
    rate_limit_sec: float
    out_dir: Path
    cache_dir: Path | None


def parse_fasta(text: str) -> List[dict]:
    """Split a FASTA payload into structured records."""
    records: List[dict] = []
    header: str | None = None
    seq_lines: list[str] = []

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if header is not None:
                records.append(
                    {
                        "header": header,
                        "accession": header.split()[0],
                        "sequence": "".join(seq_lines),
                    }
                )
            header = stripped[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(stripped)

    if header is not None:
        records.append(
            {
                "header": header,
                "accession": header.split()[0],
                "sequence": "".join(seq_lines),
            }
        )

    return records


def _cache_key(prefix: str, parts: Iterable[str]) -> str:
    return "::".join([prefix, *parts])


def _build_metadata_map(summary_payload: dict) -> dict[str, dict]:
    result = summary_payload.get("result", {}) if summary_payload else {}
    uids = result.get("uids", [])
    metadata: dict[str, dict] = {}
    for uid in uids:
        entry = result.get(uid, {})
        accession = entry.get("accessionversion") or entry.get("caption")
        if not accession:
            continue
        metadata[accession] = {
            "organism": entry.get("organism"),
            "title": entry.get("title"),
            "created": entry.get("createdate"),
            "updated": entry.get("updatedate"),
        }
    return metadata


def run_pipeline(config: RunConfig) -> dict:
    """Full end-to-end run from NCBI query to data exports."""
    client = NCBIClient(
        NCBIClientConfig(
            email=config.email,
            tool=config.tool,
            api_key=config.api_key or None,
            rate_limit_sec=config.rate_limit_sec,
        )
    )

    ids_key = _cache_key("esearch", [config.query, str(config.limit)])
    cached_ids = get_cached_text(ids_key, config.cache_dir)
    if cached_ids:
        ids = [line.strip() for line in cached_ids.splitlines() if line.strip()]
    else:
        ids = client.esearch_ids(config.query, config.limit)
        set_cached_text(ids_key, "\n".join(ids), config.cache_dir)

    if not ids:
        return {
            "query": config.query,
            "id_count": 0,
            "record_count": 0,
            "paths": {},
        }

    fetch_key = _cache_key("efetch", [",".join(ids)])
    fasta_text = get_cached_text(fetch_key, config.cache_dir)
    if not fasta_text:
        fasta_text = client.efetch_fasta(ids)
        set_cached_text(fetch_key, fasta_text, config.cache_dir)

    summary_key = _cache_key("esummary", [",".join(ids)])
    summary_text = get_cached_text(summary_key, config.cache_dir)
    if summary_text:
        summary_payload = json.loads(summary_text)
    else:
        summary_payload = client.esummary_json(ids)
        set_cached_text(summary_key, json.dumps(summary_payload), config.cache_dir)

    metadata_map = _build_metadata_map(summary_payload)

    records = parse_fasta(fasta_text)
    cleaned_records: List[dict] = []
    feature_rows: List[dict] = []
    metadata_rows: List[dict] = []

    for record in records:
        cleaned_seq = clean_sequence(record.get("sequence", ""))
        if not cleaned_seq:
            continue
        enriched_record = {
            **record,
            "sequence": cleaned_seq,
        }
        cleaned_records.append(enriched_record)
        feature = {
            "accession": enriched_record["accession"],
        }
        feature.update(compute_features(cleaned_seq))
        feature_rows.append(feature)
        meta = metadata_map.get(enriched_record["accession"], {})
        metadata_rows.append(
            {
                "accession": enriched_record["accession"],
                "header": enriched_record["header"],
                "sequence_length": len(cleaned_seq),
                "organism": meta.get("organism"),
                "title": meta.get("title"),
                "created": meta.get("created"),
                "updated": meta.get("updated"),
            }
        )

    out_dir = config.out_dir
    fasta_path = write_fasta(cleaned_records, out_dir / "sequences.fasta")
    features_path = write_csv(feature_rows, out_dir / "features.csv")
    metadata_path = write_jsonl(metadata_rows, out_dir / "metadata.jsonl")

    return {
        "query": config.query,
        "id_count": len(ids),
        "record_count": len(cleaned_records),
        "paths": {
            "fasta": str(fasta_path),
            "features": str(features_path),
            "metadata": str(metadata_path),
        },
    }
