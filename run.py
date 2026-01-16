"""Convenience runner for the GenBank to BioDataset pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path
import textwrap
import yaml

from genbank_to_biodataset.pipeline import RunConfig, run_pipeline


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Execute the GenBank-to-BioDataset pipeline without installing the CLI.",
    )
    parser.add_argument(
        "--query",
        required=True,
        help="NCBI Entrez query (ESearch syntax).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=20,
        help="Maximum number of IDs to fetch (retmax).",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="Path to configuration YAML file.",
    )
    return parser.parse_args()


def _load_config(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _resolve_path(value: str | Path | None) -> Path | None:
    if value is None:
        return None
    return Path(value).expanduser()


def main() -> int:
    args = parse_args()
    cfg = _load_config(args.config)

    email = cfg.get("email", "").strip()
    tool = cfg.get("tool", "").strip()
    if not email or not tool:
        raise ValueError("Config must define both `email` and `tool` values.")

    out_dir = _resolve_path(cfg.get("out_dir", "data/processed")) or Path("data/processed")
    cache_dir = _resolve_path(cfg.get("cache_dir"))

    run_config = RunConfig(
        query=args.query,
        limit=args.limit,
        email=email,
        tool=tool,
        api_key=(cfg.get("api_key") or "").strip() or None,
        rate_limit_sec=float(cfg.get("rate_limit_sec", 0.35)),
        out_dir=out_dir,
        cache_dir=cache_dir,
    )

    summary = run_pipeline(run_config)
    _print_summary(summary)
    return 0


def _print_summary(summary: dict) -> None:
    lines = [
        f"Query: {summary.get('query')}",
        f"IDs retrieved: {summary.get('id_count', 0)}",
        f"Records exported: {summary.get('record_count', 0)}",
    ]
    paths = summary.get("paths") or {}
    if paths:
        lines.append("Outputs:")
        for label, path in paths.items():
            lines.append(f"  - {label}: {path}")
    else:
        lines.append("No outputs generated.")
    print(textwrap.dedent("\n".join(lines)))


if __name__ == "__main__":
    raise SystemExit(main())
