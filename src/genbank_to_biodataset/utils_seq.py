"""Lightweight sequence utilities for consensus workflows."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable, List, Tuple


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    return sequence.translate(COMPLEMENT)[::-1]


def kmer_profile(sequence: str, k: int = 4) -> Counter:
    """Return a Counter of k-mers, ignoring kmers containing N."""
    counts: Counter = Counter()
    if k <= 0:
        return counts
    limit = len(sequence) - k + 1
    if limit <= 0:
        return counts
    for idx in range(limit):
        kmer = sequence[idx : idx + k]
        if "N" in kmer:
            continue
        counts[kmer] += 1
    return counts


def jaccard_similarity(a: Iterable[str] | Counter, b: Iterable[str] | Counter) -> float:
    """Compute Jaccard similarity between two iterable collections."""
    set_a = set(a if not isinstance(a, Counter) else a.keys())
    set_b = set(b if not isinstance(b, Counter) else b.keys())
    if not set_a or not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return intersection / union


def approximate_identity(seq_a: str, seq_b: str, k: int = 4) -> float:
    """Approximate similarity using Jaccard on k-mer sets."""
    profile_a = kmer_profile(seq_a, k=k)
    profile_b = kmer_profile(seq_b, k=k)
    return jaccard_similarity(profile_a, profile_b)


def load_fasta_sequences(path: Path) -> List[Tuple[str, str]]:
    """Parse a small FASTA file and return a list of (header, sequence)."""
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")
    records: List[Tuple[str, str]] = []
    header: str | None = None
    seq_lines: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = stripped[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(stripped.upper())
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records
