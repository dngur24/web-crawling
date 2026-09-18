"""
Parser: extracts structured data from the raw HTML returned by AIPpred.

AIPpred (http://www.thegleelab.org/AIPpred/) uses a Random Forest model
and returns a table with columns for peptide sequence, class (AIP / Non-AIP),
and probability score.

Expected result table columns (from live observation):
    Peptide | Class | Probability
or similar variants.
"""

from __future__ import annotations

import re
from typing import Optional

from bs4 import BeautifulSoup


class Parser:
    """
    Parses the raw HTML result page from AIPpred.

    Strategy
    --------
    1. Table-based: look for a <table> whose headers contain score/probability
       and prediction/class columns, then extract the first data row.
    2. Regex fallback: scan the plain text for recognizable patterns.
    """

    # Regex fallbacks
    _SCORE_RE = re.compile(
        r"(?:probability|score|prob)\s*[:\-]?\s*([\d.]+)", re.IGNORECASE
    )
    _LABEL_RE = re.compile(
        r"\b(Anti[- ]?Inflammatory|Non[- ]?Anti[- ]?Inflammatory|AIP|Non[- ]?AIP)\b",
        re.IGNORECASE,
    )

    def parse(self, html: str) -> tuple[Optional[float], Optional[str]]:
        """
        Extract (score, label) from the result HTML.

        Parameters
        ----------
        html : str  – raw HTML of the result page

        Returns
        -------
        (score, label) – either may be None if parsing fails.
        """
        soup = BeautifulSoup(html, "lxml")

        # ── 1. Table-based parsing ────────────────────────────────────────
        score, label = self._parse_table(soup)
        if score is not None or label is not None:
            return score, label

        # ── 2. Regex fallback on plain text ──────────────────────────────
        text = soup.get_text(separator=" ")
        score = self._regex_score(text)
        label = self._regex_label(text)
        return score, label

    # ------------------------------------------------------------------
    # Table parsing
    # ------------------------------------------------------------------

    def _parse_table(
        self, soup: BeautifulSoup
    ) -> tuple[Optional[float], Optional[str]]:
        """
        Locate a <table> whose headers contain score/class columns and
        return values from the first data row.
        """
        for table in soup.find_all("table"):
            rows = table.find_all("tr")
            if len(rows) < 2:
                continue

            header_idx, col_score, col_label = self._find_header(rows)
            if header_idx is None:
                continue

            # First data row after the header
            for row in rows[header_idx + 1 :]:
                cells = [td.get_text(strip=True) for td in row.find_all(["td", "th"])]
                if not cells:
                    continue

                score = (
                    self._safe_float(cells[col_score])
                    if col_score is not None and col_score < len(cells)
                    else None
                )
                label = (
                    self._normalize_label(cells[col_label])
                    if col_label is not None and col_label < len(cells) and cells[col_label]
                    else None
                )
                return score, label

        return None, None

    def _find_header(
        self, rows: list
    ) -> tuple[Optional[int], Optional[int], Optional[int]]:
        """
        Search rows for a header containing score/probability and class/prediction columns.

        Actual AIPpred columns: S.NO | FASTA ID | AIP or Non-AIP | Prob

        Returns (header_row_idx, col_score, col_label).
        """
        for i, row in enumerate(rows):
            cells = [td.get_text(strip=True).lower() for td in row.find_all(["td", "th"])]

            has_score = any(
                kw in c for c in cells for kw in ("prob", "score")
            )
            has_pred = any(
                kw in c for c in cells for kw in ("aip", "class", "predict", "result")
            )

            if not (has_score or has_pred):
                continue

            col_score = next(
                (j for j, c in enumerate(cells) if "prob" in c or "score" in c), None
            )
            # "AIP or Non-AIP" column → contains "aip"
            col_label = next(
                (
                    j
                    for j, c in enumerate(cells)
                    if "aip" in c or "class" in c or "predict" in c or "result" in c
                ),
                None,
            )
            return i, col_score, col_label

        return None, None, None

    # ------------------------------------------------------------------
    # Regex helpers
    # ------------------------------------------------------------------

    def _regex_score(self, text: str) -> Optional[float]:
        m = self._SCORE_RE.search(text)
        return self._safe_float(m.group(1)) if m else None

    def _regex_label(self, text: str) -> Optional[str]:
        m = self._LABEL_RE.search(text)
        return self._normalize_label(m.group(0)) if m else None

    # ------------------------------------------------------------------
    # Static helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_float(value: str) -> Optional[float]:
        try:
            return float(value)
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _normalize_label(label: str) -> str:
        """
        Standardize label variants from the server.

        The server returns 'AIP' or 'Non-AIP'. Both are mapped to
        human-readable forms for consistency with the rest of the project.
        """
        lower = label.lower().replace("-", "").replace(" ", "")
        if "non" in lower:
            return "Non Anti-Inflammatory"
        # Catches 'AIP', 'Anti-Inflammatory', 'Anti Inflammatory', etc.
        return "Anti-Inflammatory"
