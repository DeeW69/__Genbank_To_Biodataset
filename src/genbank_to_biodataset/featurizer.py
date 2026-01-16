"""Sequence cleaning and feature computation."""

from __future__ import annotations

ALLOWED_BASES = {"A", "C", "G", "T", "N"}


def clean_sequence(sequence: str) -> str:
    """Uppercase and keep only canonical bases."""
    upper = sequence.upper()
    return "".join(ch for ch in upper if ch in ALLOWED_BASES)


def compute_features(sequence: str) -> dict:
    """Return simple descriptive statistics for a nucleotide sequence."""
    length = len(sequence)
    counts = {base: sequence.count(base) for base in ALLOWED_BASES}
    gc = counts["G"] + counts["C"]
    gc_percent = (gc / length * 100) if length else 0.0
    return {
        "length": length,
        "gc_percent": round(gc_percent, 4),
        "count_A": counts["A"],
        "count_C": counts["C"],
        "count_G": counts["G"],
        "count_T": counts["T"],
        "count_N": counts["N"],
    }
