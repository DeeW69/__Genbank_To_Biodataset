"""Simple text cache utilities using SHA256 keys."""

from __future__ import annotations

import hashlib
from pathlib import Path


def _key_to_path(cache_dir: Path | None, key: str) -> Path | None:
    if cache_dir is None:
        return None
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return cache_dir / f"{digest}.txt"


def get_cached_text(key: str, cache_dir: Path | None) -> str | None:
    path = _key_to_path(cache_dir, key)
    if path is None or not path.exists():
        return None
    return path.read_text(encoding="utf-8")


def set_cached_text(key: str, value: str, cache_dir: Path | None) -> None:
    path = _key_to_path(cache_dir, key)
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8")
