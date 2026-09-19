"""
Parser: extracts structured data from the raw HTML returned by ACVPICPred.

Target: https://i.uestc.edu.cn/acvpICPred/main/Main.php

The result page is expected to contain a table with predicted activity values
and classification labels for the submitted peptide.
"""

from __future__ import annotations

import re
from typing import Optional

from bs4 import BeautifulSoup


class Parser:
    """
    Parses the raw HTML result page from ACVPICPred.

    The server returns a result page containing:
    - A predicted numerical value (IC50 / EC50 / IC90 in the selected unit)
    - A classification label (e.g. "Active" / "Inactive")

    Strategy
    --------
    1. Table-based parsing: find a <table> whose headers contain
       value/score/activity and label/prediction columns.
    2. Regex fallback: scan plain text for numeric values and label keywords.
    """

    # ── Regex fallbacks ───────────────────────────────────────────────────────
    _VALUE_RE = re.compile(
        r"(?:predicted|IC50|EC50|IC90|activity|value)\s*[:\-]?\s*([-\d.]+)",
        re.IGNORECASE,
    )
    _BARE_FLOAT_RE = re.compile(r"([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)")

    _LABEL_RE = re.compile(
        r"\b(Active|Inactive|Anti-?coronavirus|Non[- ]?Anti-?coronavirus)\b",
        re.IGNORECASE,
    )

    # Column header keywords
    _VALUE_KEYWORDS  = ("value", "score", "activity", "predicted", "ic50", "ec50", "ic90")
    _LABEL_KEYWORDS  = ("label", "predict", "result", "class", "activity class")

    # ─────────────────────────────────────────────────────────────────────────

    def parse(self, html: str) -> tuple[Optional[float], Optional[str]]:
        """
        Extract ``(predicted_value, label)`` from the result HTML.

        Parameters
        ----------
        html : str
            Raw HTML of the ACVPICPred result page.

        Returns
        -------
        (predicted_value, label)
            Either element may be ``None`` if parsing fails.
        """
        soup = BeautifulSoup(html, "lxml")

        # ── 1. Table-based parsing ────────────────────────────────────────────
        value, label = self._parse_table(soup)
        if value is not None or label is not None:
            return value, label

        # ── 2. Regex fallback on visible text ────────────────────────────────
        text = soup.get_text(separator=" ")
        value = self._regex_value(text)
        label = self._regex_label(text)
        return value, label

    # ------------------------------------------------------------------
    # Table parsing
    # ------------------------------------------------------------------

    def _parse_table(
        self, soup: BeautifulSoup
    ) -> tuple[Optional[float], Optional[str]]:
        """
        Find a result table and extract the predicted value and label.

        Returns (predicted_value, label).
        """
        for table in soup.find_all("table"):
            rows = table.find_all("tr")
            if len(rows) < 2:
                continue

            header_idx, col_value, col_label = self._find_header(rows)
            if header_idx is None:
                continue

            # Take the first data row after the header
            for row in rows[header_idx + 1:]:
                cells = [td.get_text(strip=True) for td in row.find_all(["td", "th"])]
                if not cells:
                    continue

                value = (
                    self._safe_float(cells[col_value])
                    if col_value is not None and col_value < len(cells)
                    else None
                )
                label = (
                    self._normalize_label(cells[col_label])
                    if col_label is not None and col_label < len(cells) and cells[col_label]
                    else None
                )
                return value, label

        return None, None

    def _find_header(
        self, rows: list
    ) -> tuple[Optional[int], Optional[int], Optional[int]]:
        """
        Search rows for a header containing value/label columns.

        Returns (header_row_idx, col_value, col_label).
        """
        for i, row in enumerate(rows):
            cells_lower = [
                td.get_text(strip=True).lower()
                for td in row.find_all(["td", "th"])
            ]

            has_value = any(
                any(kw in c for kw in self._VALUE_KEYWORDS) for c in cells_lower
            )
            has_label = any(
                any(kw in c for kw in self._LABEL_KEYWORDS) for c in cells_lower
            )

            if not (has_value or has_label):
                continue

            col_value = next(
                (j for j, c in enumerate(cells_lower)
                 if any(kw in c for kw in self._VALUE_KEYWORDS)),
                None,
            )
            col_label = next(
                (j for j, c in enumerate(cells_lower)
                 if any(kw in c for kw in self._LABEL_KEYWORDS)),
                None,
            )
            return i, col_value, col_label

        return None, None, None

    # ------------------------------------------------------------------
    # Regex helpers
    # ------------------------------------------------------------------

    def _regex_value(self, text: str) -> Optional[float]:
        """Try to extract a numeric predicted value from plain text."""
        m = self._VALUE_RE.search(text)
        if m:
            return self._safe_float(m.group(1))
        # Broader fallback: grab first float/int that looks like a result
        m2 = self._BARE_FLOAT_RE.search(text)
        return self._safe_float(m2.group(1)) if m2 else None

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
        """Standardize label text from the server."""
        lower = label.lower()
        if "non" in lower or "inactive" in lower:
            return "Inactive"
        return "Active"
