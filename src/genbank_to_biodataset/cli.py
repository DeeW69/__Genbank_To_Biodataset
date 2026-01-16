"""Command line interface for genbank-to-biodataset."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import textwrap
import yaml

from .pipeline import RunConfig, run_pipeline
from .consensus import ConsensusConfig, run_consensus


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genbank2bio",
        description="Build ML-ready datasets directly from NCBI GenBank.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run",
        help="Execute the full download/processing pipeline.",
    )
    run_parser.add_argument(
        "--query",
        required=True,
        help="NCBI query string (ESearch syntax).",
    )
    run_parser.add_argument(
        "--limit",
        type=int,
        default=50,
        help="Maximum number of IDs to fetch (retmax).",
    )
    run_parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="Path to the YAML configuration file.",
    )
    consensus_parser = subparsers.add_parser(
        "consensus",
        help="Construire des consensus experimentaux (COI/COX1) a partir de fragments.",
    )
    consensus_parser.add_argument(
        "--query",
        required=True,
        help="Requete Entrez pour la construction des consensus.",
    )
    consensus_parser.add_argument(
        "--limit",
        type=int,
        default=200,
        help="Nombre maximal d'IDs a recuperer (retmax).",
    )
    consensus_parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="Chemin vers le fichier de configuration YAML.",
    )
    consensus_parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/processed_consensus"),
        help="Dossier de sortie pour les consensus (defaut: data/processed_consensus).",
    )
    consensus_parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("data/cache"),
        help="Dossier cache local (defaut: data/cache).",
    )
    consensus_parser.add_argument(
        "--min-len",
        type=int,
        default=200,
        help="Longueur minimale des fragments conserves.",
    )
    consensus_parser.add_argument(
        "--max-n-frac",
        type=float,
        default=0.05,
        help="Fraction maximale de N autorisee (ex: 0.05 = 5%%).",
    )
    consensus_parser.add_argument(
        "--min-overlap",
        type=int,
        default=80,
        help="Chevauchement minimal pour fusionner deux fragments.",
    )
    consensus_parser.add_argument(
        "--min-identity",
        type=float,
        default=0.97,
        help="Identite minimale pour rattacher une sequence a un cluster.",
    )
    consensus_parser.add_argument(
        "--max-seqs-per-cluster",
        type=int,
        default=200,
        help="Nombre maximum de sequences par cluster (pour eviter la deroute).",
    )
    consensus_parser.add_argument(
        "--reference-fasta",
        type=Path,
        help="FASTA de reference pour orienter les sequences (optionnel).",
    )
    consensus_parser.add_argument(
        "--use-mafft",
        action="store_true",
        help="Utiliser MAFFT si disponible pour construire le consensus.",
    )
    consensus_parser.add_argument(
        "--use-vsearch",
        action="store_true",
        help="Utiliser vsearch pour le clustering si installe.",
    )
    return parser


def _load_config(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _resolve_path(value: str | Path | None) -> Path | None:
    if value is None:
        return None
    return Path(value).expanduser()


def _run_command(args: argparse.Namespace) -> int:
    raw_config = _load_config(args.config)
    email = raw_config.get("email", "").strip()
    tool = raw_config.get("tool", "").strip()
    if not email:
        raise ValueError("`email` must be set in the configuration file.")
    if not tool:
        raise ValueError("`tool` must be set in the configuration file.")

    out_dir = _resolve_path(raw_config.get("out_dir", "data/processed")) or Path("data/processed")
    cache_dir = _resolve_path(raw_config.get("cache_dir"))

    run_config = RunConfig(
        query=args.query,
        limit=args.limit,
        email=email,
        tool=tool,
        api_key=(raw_config.get("api_key") or "").strip() or None,
        rate_limit_sec=float(raw_config.get("rate_limit_sec", 0.35)),
        out_dir=out_dir,
        cache_dir=cache_dir,
    )

    summary = run_pipeline(run_config)
    _print_summary(summary)
    return 0


def _run_consensus_command(args: argparse.Namespace) -> int:
    raw_config = _load_config(args.config)
    email = raw_config.get("email", "").strip()
    tool = raw_config.get("tool", "").strip()
    if not email:
        raise ValueError("`email` doit etre defini dans la configuration.")
    if not tool:
        raise ValueError("`tool` doit etre defini dans la configuration.")

    run_config = ConsensusConfig(
        query=args.query,
        limit=args.limit,
        email=email,
        tool=tool,
        api_key=(raw_config.get("api_key") or "").strip() or None,
        rate_limit_sec=float(raw_config.get("rate_limit_sec", 0.35)),
        out_dir=_resolve_path(args.out_dir) or Path("data/processed_consensus"),
        cache_dir=_resolve_path(args.cache_dir),
        min_len=args.min_len,
        max_n_frac=args.max_n_frac,
        min_overlap=args.min_overlap,
        min_identity=args.min_identity,
        max_seqs_per_cluster=args.max_seqs_per_cluster,
        reference_fasta=_resolve_path(args.reference_fasta),
        use_mafft=bool(args.use_mafft),
        use_vsearch=bool(args.use_vsearch),
    )

    summary = run_consensus(run_config)
    _print_consensus_summary(summary)
    return 0


def _print_summary(summary: dict) -> None:
    lines = [
        f"Query: {summary.get('query')}",
        f"IDs retrieved: {summary.get('id_count', 0)}",
        f"Records exported: {summary.get('record_count', 0)}",
    ]

    if summary.get("paths"):
        lines.append("Outputs:")
        for label, path in summary["paths"].items():
            lines.append(f"  - {label}: {path}")
    else:
        lines.append("No records were exported.")

    message = "\n".join(lines)
    print(textwrap.dedent(message))


def _print_consensus_summary(summary: dict) -> None:
    lines = [
        f"Query: {summary.get('query')}",
        f"IDs recuperes: {summary.get('id_count', 0)}",
        f"Sequences filtrees: {summary.get('filtered_count', 0)}",
        f"Clusters: {summary.get('clusters', 0)}",
        f"Consensus experimentaux: {summary.get('consensus', 0)}",
    ]
    paths = summary.get("paths") or {}
    if paths:
        lines.append("Fichiers generes:")
        for label, path in paths.items():
            lines.append(f"  - {label}: {path}")
    print(textwrap.dedent("\n".join(lines)))


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "run":
        return _run_command(args)
    if args.command == "consensus":
        return _run_consensus_command(args)
    parser.print_help()
    return 1


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
