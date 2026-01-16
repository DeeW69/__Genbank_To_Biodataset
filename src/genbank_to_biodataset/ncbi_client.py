"""NCBI E-utilities client built on top of requests."""

from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Iterable, List

import requests

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(slots=True)
class NCBIClientConfig:
    email: str
    tool: str
    api_key: str | None = None
    rate_limit_sec: float = 0.35


class NCBIClient:
    """Minimal helper for ESearch/EFetch/ESummary calls with rate limiting."""

    def __init__(self, config: NCBIClientConfig) -> None:
        if not config.email:
            raise ValueError("NCBI email must be provided.")
        if not config.tool:
            raise ValueError("NCBI tool value must be provided.")
        self._config = config
        self._session = requests.Session()
        self._last_call = 0.0

    def _throttle(self) -> None:
        now = time.monotonic()
        delta = now - self._last_call
        wait = max(0.0, self._config.rate_limit_sec - delta)
        if wait:
            time.sleep(wait)

    def _request(self, endpoint: str, params: dict) -> requests.Response:
        self._throttle()
        query = {
            **params,
            "email": self._config.email,
            "tool": self._config.tool,
        }
        if self._config.api_key:
            query["api_key"] = self._config.api_key

        url = f"{BASE_URL.rstrip('/')}/{endpoint.lstrip('/')}"
        response = self._session.get(url, params=query, timeout=60)
        response.raise_for_status()
        self._last_call = time.monotonic()
        return response

    def esearch_ids(self, query: str, limit: int) -> List[str]:
        """Run an ESearch query against nuccore and return ID list."""
        response = self._request(
            "esearch.fcgi",
            {
                "db": "nuccore",
                "term": query,
                "retmode": "json",
                "retmax": limit,
            },
        )
        payload = response.json()
        id_list = payload.get("esearchresult", {}).get("idlist", [])
        return list(id_list)

    def efetch_fasta(self, ids: Iterable[str]) -> str:
        """Fetch FASTA records for the provided IDs."""
        id_param = ",".join([i for i in ids if i])
        if not id_param:
            return ""
        response = self._request(
            "efetch.fcgi",
            {
                "db": "nuccore",
                "id": id_param,
                "rettype": "fasta",
                "retmode": "text",
            },
        )
        return response.text

    def esummary_json(self, ids: Iterable[str]) -> dict:
        """Return ESummary JSON for the provided IDs."""
        id_param = ",".join([i for i in ids if i])
        if not id_param:
            return {}
        response = self._request(
            "esummary.fcgi",
            {
                "db": "nuccore",
                "id": id_param,
                "retmode": "json",
            },
        )
        return response.json()
