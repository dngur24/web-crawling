"""
Parser: extracts structured data from the raw HTML returned by AagingPEPred.

The server (process_pep_pred.php) returns an HTML page containing a results
table.  Typical column headers observed:
    Peptide ID | Sequence | Score | Prediction

This parser handles:
 1. Table-based extraction (primary)
 2. Regex fallback on plain text
"""

from __future__ import annotations

import re
from typing import Optional

from bs4 import BeautifulSoup


class Parser:
    """
    Parses the raw HTML result page from AagingPEPred.

    Returns (score, label) for the first (or only) result row found.
    """

    # ── Regex fallbacks ──────────────────────────────────────────────────
    _SCORE_RE = re.compile(r"score\s*[:\-]?\s*([-\d.]+)", re.IGNORECASE)
    _LABEL_RE = re.compile(
        r"(Anti[- ]?[Aa]ging|Non[- ]?Anti[- ]?[Aa]ging|Non[- ]?[Aa]ging)",
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
        Locate a table with score/prediction headers and extract values
        from the first data row.
        """
        for table in soup.find_all("table"):
            rows = table.find_all("tr")
            if len(rows) < 2:
                continue

            header_idx, col_score, col_label = self._find_header(rows)
            if header_idx is None:
                continue

            # Return values from the first data row
            for row in rows[header_idx + 1:]:
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
        Search rows for a header containing score/prediction columns.

        AagingPEPred uses: 'XGB_Score' for score, 'Class' for label.
        Fallback keywords also handled for robustness.

        Returns (header_row_idx, col_score, col_label).
        """
        for i, row in enumerate(rows):
            cells = [td.get_text(strip=True).lower() for td in row.find_all(["td", "th"])]

            # AagingPEPred-specific headers: XGB_Score, Class
            has_score = any("score" in c for c in cells)
            has_pred = any(
                c in ("class",) or "predict" in c or "result" in c or "label" in c
                for c in cells
            )
            if not (has_score or has_pred):
                continue

            col_score = next(
                (j for j, c in enumerate(cells) if "score" in c), None
            )
            col_label = next(
                (j for j, c in enumerate(cells)
                 if c == "class" or "predict" in c or "result" in c or "label" in c),
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

        Known server values: 'Anti-aging', 'Non Anti-aging'
        """
        lower = label.lower()
        if "non" in lower:
            return "Non Anti-Aging"
        if "anti" in lower or "aging" in lower:
            return "Anti-Aging"
        return label  # pass through unknown values unchanged
