"""Experimental consensus reconstruction pipeline."""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

from .cache import get_cached_text, set_cached_text
from .exporter import write_csv, write_fasta, write_jsonl
from .featurizer import clean_sequence, compute_features
from .ncbi_client import NCBIClient, NCBIClientConfig
from .pipeline import parse_fasta
from .utils_seq import (
    approximate_identity,
    load_fasta_sequences,
    reverse_complement,
)

LOGGER = logging.getLogger(__name__)


@dataclass(slots=True)
class ConsensusConfig:
    query: str
    limit: int
    email: str
    tool: str
    api_key: str | None
    rate_limit_sec: float
    out_dir: Path
    cache_dir: Path | None
    min_len: int
    max_n_frac: float
    min_overlap: int
    min_identity: float
    max_seqs_per_cluster: int
    reference_fasta: Path | None
    use_mafft: bool
    use_vsearch: bool


def run_consensus(config: ConsensusConfig) -> dict:
    """Main orchestration for consensus workflow."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    client = NCBIClient(
        NCBIClientConfig(
            email=config.email,
            tool=config.tool,
            api_key=config.api_key or None,
            rate_limit_sec=config.rate_limit_sec,
        )
    )

    ids = _load_ids(client, config)
    fasta_text = _load_fasta(client, config, ids)
    metadata_map = _load_metadata(client, config, ids)

    raw_records = parse_fasta(fasta_text)
    LOGGER.info("Sequences retrieved: %s", len(raw_records))

    reference_seq = _load_reference_sequence(config.reference_fasta)
    filtered_records = _filter_records(raw_records, config, metadata_map, reference_seq)
    LOGGER.info("Sequences retained after QC: %s", len(filtered_records))

    if not filtered_records:
        return {
            "query": config.query,
            "id_count": len(ids),
            "filtered_count": 0,
            "clusters": 0,
            "consensus": 0,
            "paths": {},
        }

    clusters = _cluster_records(filtered_records, config, reference_seq)
    LOGGER.info("Clusters created: %s", len(clusters))

    consensus_entries, cluster_stats = _build_consensus_entries(clusters, config)
    LOGGER.info("Consensus sequences generated: %s", len(consensus_entries))

    cluster_json_rows = _clusters_to_json(clusters, cluster_stats)
    filtered_metadata_rows, filtered_feature_rows = _build_filtered_outputs(
        clusters, metadata_map
    )

    outputs = {
        "consensus_fasta": write_fasta(
            consensus_entries, config.out_dir / "consensus.fasta"
        ),
        "consensus_stats": write_csv(
            cluster_stats, config.out_dir / "consensus_stats.csv"
        ),
        "clusters": write_jsonl(
            cluster_json_rows, config.out_dir / "clusters.jsonl"
        ),
        "filtered_metadata": write_jsonl(
            filtered_metadata_rows, config.out_dir / "filtered_metadata.jsonl"
        ),
        "filtered_features": write_csv(
            filtered_feature_rows, config.out_dir / "filtered_features.csv"
        ),
    }

    return {
        "query": config.query,
        "id_count": len(ids),
        "filtered_count": sum(len(cluster["records"]) for cluster in clusters),
        "clusters": len(clusters),
        "consensus": len(consensus_entries),
        "paths": {key: str(path) for key, path in outputs.items()},
    }


def _load_ids(client: NCBIClient, config: ConsensusConfig) -> List[str]:
    cache_key = _cache_key("esearch", [config.query, str(config.limit)])
    cached = get_cached_text(cache_key, config.cache_dir)
    if cached:
        return [line.strip() for line in cached.splitlines() if line.strip()]
    ids = client.esearch_ids(config.query, config.limit)
    set_cached_text(cache_key, "\n".join(ids), config.cache_dir)
    return ids


def _load_fasta(client: NCBIClient, config: ConsensusConfig, ids: Iterable[str]) -> str:
    cache_key = _cache_key("efetch", [",".join(ids)])
    cached = get_cached_text(cache_key, config.cache_dir)
    if cached:
        return cached
    fasta = client.efetch_fasta(ids)
    set_cached_text(cache_key, fasta, config.cache_dir)
    return fasta


def _load_metadata(
    client: NCBIClient, config: ConsensusConfig, ids: Iterable[str]
) -> dict:
    cache_key = _cache_key("esummary", [",".join(ids)])
    cached = get_cached_text(cache_key, config.cache_dir)
    if cached:
        payload = json.loads(cached)
    else:
        payload = client.esummary_json(ids)
        set_cached_text(cache_key, json.dumps(payload), config.cache_dir)
    result = payload.get("result", {}) if payload else {}
    metadata: dict[str, dict] = {}
    for uid in result.get("uids", []):
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


def _load_reference_sequence(path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        records = load_fasta_sequences(path)
    except FileNotFoundError:
        LOGGER.warning("Reference FASTA introuvable: %s", path)
        return None
    if not records:
        LOGGER.warning("Reference FASTA vide: %s", path)
        return None
    return records[0][1]


def _filter_records(
    records: List[dict],
    config: ConsensusConfig,
    metadata_map: dict,
    reference_seq: str | None,
) -> List[dict]:
    filtered: List[dict] = []
    for record in records:
        seq = clean_sequence(record.get("sequence", ""))
        if not seq:
            continue
        if len(seq) < config.min_len:
            continue
        n_frac = seq.count("N") / len(seq)
        if n_frac > config.max_n_frac:
            continue
        if reference_seq:
            seq = _orient_to_reference(seq, reference_seq)
        filtered.append(
            {
                "accession": record["accession"],
                "header": record["header"],
                "sequence": seq,
                "metadata": metadata_map.get(record["accession"], {}),
            }
        )
    return filtered


def _orient_to_reference(sequence: str, reference_seq: str) -> str:
    direct = approximate_identity(sequence, reference_seq)
    flipped = approximate_identity(reverse_complement(sequence), reference_seq)
    return reverse_complement(sequence) if flipped > direct else sequence


def _cluster_records(
    records: List[dict],
    config: ConsensusConfig,
    reference_seq: str | None,
) -> List[dict]:
    if config.use_vsearch and _binary_available("vsearch"):
        try:
            return _cluster_with_vsearch(records, config)
        except Exception as exc:  # pragma: no cover - defensive
            LOGGER.warning("vsearch a echoue (%s), repli sur clustering Python.", exc)
    elif config.use_vsearch:
        LOGGER.warning("vsearch non detecte, clustering Python utilise.")
    return _cluster_greedy(records, config, reference_seq)


def _cluster_greedy(
    records: List[dict],
    config: ConsensusConfig,
    reference_seq: str | None,
) -> List[dict]:
    remaining = [dict(record) for record in records]
    clusters: List[dict] = []
    cluster_idx = 1
    while remaining:
        seed = remaining.pop(0)
        cluster_records = [seed]
        seed_seq = seed["sequence"]
        assigned_indices: List[int] = []
        for idx, candidate in enumerate(remaining):
            seq = candidate["sequence"]
            oriented_seq, identity = _orient_to_seed(seq, seed_seq)
            if identity >= config.min_identity:
                candidate["sequence"] = oriented_seq
                cluster_records.append(candidate)
                assigned_indices.append(idx)
                if len(cluster_records) >= config.max_seqs_per_cluster:
                    break
        for offset, rem_idx in enumerate(assigned_indices):
            remaining.pop(rem_idx - offset)
        clusters.append(
            {
                "cluster_id": f"cluster_{cluster_idx}",
                "records": cluster_records,
            }
        )
        cluster_idx += 1
    return clusters


def _cluster_with_vsearch(records: List[dict], config: ConsensusConfig) -> List[dict]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        fasta_path = tmp_path / "input.fasta"
        uc_path = tmp_path / "clusters.uc"
        with fasta_path.open("w", encoding="utf-8") as handle:
            for idx, record in enumerate(records):
                handle.write(f">seq{idx}|{record['accession']}\n")
                handle.write(f"{record['sequence']}\n")

        cmd = [
            "vsearch",
            "--cluster_fast",
            str(fasta_path),
            "--id",
            f"{config.min_identity:.2f}",
            "--strand",
            "both",
            "--uc",
            str(uc_path),
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        id_map = {f"seq{idx}|{record['accession']}": record for idx, record in enumerate(records)}
        clusters: dict[str, List[dict]] = {}
        with uc_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                cluster_num = parts[1]
                strand = parts[4]
                label = parts[8]
                source = dict(id_map.get(label, {}))
                if not source:
                    continue
                if strand == "-":
                    source["sequence"] = reverse_complement(source["sequence"])
                clusters.setdefault(cluster_num, []).append(source)

        cluster_list: List[dict] = []
        cluster_idx = 1
        for _, seqs in clusters.items():
            if not seqs:
                continue
            for chunk in _chunk(seqs, config.max_seqs_per_cluster):
                cluster_list.append(
                    {
                        "cluster_id": f"cluster_{cluster_idx}",
                        "records": chunk,
                    }
                )
                cluster_idx += 1
        return cluster_list


def _chunk(data: List[dict], size: int) -> Iterable[List[dict]]:
    for idx in range(0, len(data), size):
        yield data[idx : idx + size]


def _orient_to_seed(sequence: str, seed_seq: str) -> tuple[str, float]:
    identity_direct = approximate_identity(sequence, seed_seq)
    rev = reverse_complement(sequence)
    identity_rev = approximate_identity(rev, seed_seq)
    if identity_rev > identity_direct:
        return rev, identity_rev
    return sequence, identity_direct


def _build_consensus_entries(
    clusters: List[dict],
    config: ConsensusConfig,
) -> tuple[List[dict], List[dict]]:
    consensus_records: List[dict] = []
    stats_rows: List[dict] = []

    mafft_ok = config.use_mafft and _binary_available("mafft")
    if config.use_mafft and not mafft_ok:
        LOGGER.warning("MAFFT non detecte, fusion overlap Python utilisee.")

    for cluster in clusters:
        sequences = [record["sequence"] for record in cluster["records"]]
        if not sequences:
            continue
        if mafft_ok and len(sequences) > 1:
            consensus_seq, support = _consensus_with_mafft(sequences)
        else:
            consensus_seq, support = _consensus_overlap(sequences, config.min_overlap)
        features = compute_features(consensus_seq)
        length = features["length"]
        n_frac = (features["count_N"] / length) if length else 0.0
        positive_support = [value for value in support if value > 0]
        mean_support = (
            sum(positive_support) / len(positive_support) if positive_support else 0.0
        )
        header = (
            f"{cluster['cluster_id']}|size={len(sequences)}|len={length}|"
            f"N%={round(n_frac * 100, 2)}|support={round(mean_support, 3)}|experimental"
        )
        consensus_records.append(
            {
                "header": header,
                "sequence": consensus_seq,
            }
        )
        stats_rows.append(
            {
                "cluster_id": cluster["cluster_id"],
                "n_seqs": len(sequences),
                "length": length,
                "gc_percent": features["gc_percent"],
                "n_frac": round(n_frac, 5),
                "mean_support": round(mean_support, 5),
            }
        )
    return consensus_records, stats_rows


def _clusters_to_json(clusters: List[dict], stats_rows: List[dict]) -> List[dict]:
    stats_map = {row["cluster_id"]: row for row in stats_rows}
    json_rows: List[dict] = []
    for cluster in clusters:
        seq_lengths = [len(rec["sequence"]) for rec in cluster["records"]]
        json_rows.append(
            {
                "cluster_id": cluster["cluster_id"],
                "size": len(cluster["records"]),
                "accessions": [rec["accession"] for rec in cluster["records"]],
                "length_min": min(seq_lengths) if seq_lengths else 0,
                "length_max": max(seq_lengths) if seq_lengths else 0,
                "mean_support": stats_map.get(cluster["cluster_id"], {}).get(
                    "mean_support"
                ),
            }
        )
    return json_rows


def _build_filtered_outputs(
    clusters: List[dict],
    metadata_map: dict,
) -> tuple[List[dict], List[dict]]:
    metadata_rows: List[dict] = []
    feature_rows: List[dict] = []
    for cluster in clusters:
        for record in cluster["records"]:
            accession = record["accession"]
            metadata = metadata_map.get(accession, {})
            metadata_rows.append(
                {
                    "cluster_id": cluster["cluster_id"],
                    "accession": accession,
                    "header": record["header"],
                    "sequence_length": len(record["sequence"]),
                    "organism": metadata.get("organism"),
                    "title": metadata.get("title"),
                    "created": metadata.get("created"),
                    "updated": metadata.get("updated"),
                }
            )
            features = compute_features(record["sequence"])
            feature_rows.append(
                {
                    "cluster_id": cluster["cluster_id"],
                    "accession": accession,
                    **features,
                }
            )
    return metadata_rows, feature_rows


def _consensus_overlap(
    sequences: List[str], min_overlap: int
) -> tuple[str, List[float]]:
    counts: List[dict[str, int]] = []
    if not sequences:
        return "", []
    for base in sequences[0]:
        counts.append({base: 1})

    for seq in sequences[1:]:
        consensus_seq, _ = _counts_to_sequence(counts)
        placement = _best_overlap(consensus_seq, seq, min_overlap)
        if placement is None:
            start_idx = len(consensus_seq) + min_overlap
        else:
            direction, overlap = placement
            if direction == "append":
                start_idx = len(consensus_seq) - overlap
            else:
                start_idx = -(len(seq) - overlap)
        _apply_sequence(counts, seq, start_idx)

    return _counts_to_sequence(counts)


def _consensus_with_mafft(sequences: List[str]) -> tuple[str, List[float]]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        input_fasta = tmp_path / "input.fasta"
        output_fasta = tmp_path / "aligned.fasta"
        with input_fasta.open("w", encoding="utf-8") as handle:
            for idx, seq in enumerate(sequences):
                handle.write(f">seq{idx}\n{seq}\n")
        cmd = ["mafft", "--auto", str(input_fasta)]
        with output_fasta.open("w", encoding="utf-8") as out_handle:
            subprocess.run(
                cmd,
                check=True,
                stdout=out_handle,
                stderr=subprocess.PIPE,
            )
        aligned = load_fasta_sequences(output_fasta)

    columns = list(zip(*[seq for _, seq in aligned]))
    consensus_chars: List[str] = []
    support: List[float] = []
    for column in columns:
        bases = [base for base in column if base != "-"]
        if not bases:
            consensus_chars.append("N")
            support.append(0.0)
            continue
        counts = {}
        for base in bases:
            counts[base] = counts.get(base, 0) + 1
        base, count = max(counts.items(), key=lambda item: item[1])
        winners = [b for b, v in counts.items() if v == count]
        if len(winners) > 1:
            consensus_chars.append("N")
        else:
            consensus_chars.append(base)
        support.append(round(count / len(bases), 5))
    consensus_seq = "".join(consensus_chars).replace("-", "")
    return consensus_seq, support


def _best_overlap(consensus: str, sequence: str, min_overlap: int):
    best = None
    best_score = -1.0
    max_overlap = min(len(consensus), len(sequence))
    for overlap in range(max_overlap, min_overlap - 1, -1):
        suffix = consensus[-overlap:]
        prefix = sequence[:overlap]
        score = _overlap_score(suffix, prefix)
        if score > best_score:
            best_score = score
            best = ("append", overlap)
        prefix_consensus = consensus[:overlap]
        suffix_seq = sequence[-overlap:]
        score_rev = _overlap_score(prefix_consensus, suffix_seq)
        if score_rev > best_score:
            best_score = score_rev
            best = ("prepend", overlap)
    return best


def _overlap_score(a: str, b: str) -> float:
    if not a or not b or len(a) != len(b):
        return 0.0
    matches = 0
    effective = 0
    for base_a, base_b in zip(a, b):
        if base_a == "N" and base_b == "N":
            continue
        effective += 1
        if base_a == base_b:
            matches += 1
    if effective == 0:
        return 0.0
    return matches / effective


def _apply_sequence(counts: List[dict[str, int]], sequence: str, start_idx: int) -> None:
    if start_idx < 0:
        counts[:0] = [{} for _ in range(-start_idx)]
        start_idx = 0
    if start_idx > len(counts):
        counts.extend({} for _ in range(start_idx - len(counts)))
    required = start_idx + len(sequence)
    if required > len(counts):
        counts.extend({} for _ in range(required - len(counts)))
    for idx, base in enumerate(sequence):
        slot = counts[start_idx + idx]
        slot[base] = slot.get(base, 0) + 1


def _counts_to_sequence(counts: List[dict[str, int]]) -> tuple[str, List[float]]:
    sequence_chars: List[str] = []
    support: List[float] = []
    for slot in counts:
        if not slot:
            sequence_chars.append("N")
            support.append(0.0)
            continue
        total = sum(slot.values())
        base, count = max(slot.items(), key=lambda item: item[1])
        winners = [b for b, v in slot.items() if v == count]
        if len(winners) > 1:
            sequence_chars.append("N")
        else:
            sequence_chars.append(base)
        support.append(round(count / total, 5))
    return "".join(sequence_chars), support


def _cache_key(prefix: str, parts: Iterable[str]) -> str:
    return "::".join([prefix, *parts])


def _binary_available(name: str) -> bool:
    return shutil.which(name) is not None
