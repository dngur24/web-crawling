"""
Parser: extracts structured data from the raw HTML returned by AntiAngioPred.

Result page table structure (from live observation):
    Original Peptide | Peptide Sequence | Mutation Position | SVM score | Prediction | ...
The first data row (Mutation Position = "No Mutation") is the prediction for the
original peptide without any mutation.
"""

from __future__ import annotations
import re
from typing import Optional
from bs4 import BeautifulSoup


class Parser:
    """
    Parses the raw HTML result page from AntiAngioPred.

    The server returns a page with a table whose headers include
    'SVM score' and 'Prediction'. The first data row (No Mutation)
    is the original peptide prediction.
    """

    # Regex fallbacks
    _SCORE_RE = re.compile(r"SVM\s+score\s*[:\-]?\s*([-\d.]+)", re.IGNORECASE)
    _SCORE_RE2 = re.compile(r"score\s*[:\-]?\s*([-\d.]+)", re.IGNORECASE)
    _LABEL_RE = re.compile(
        r"(Anti-angiogenic|Non[- ]?Anti-angiogenic|Angiogenic|Non[- ]?Angiogenic)",
        re.IGNORECASE,
    )

    def parse(self, html: str) -> tuple[Optional[float], Optional[str]]:
        """
        Extract (score, label) from the result HTML.

        Returns the prediction for the original (non-mutated) peptide.

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
        Locate a table with 'SVM score' / 'Prediction' headers and
        return values from the 'No Mutation' (original peptide) row.
        """
        for table in soup.find_all("table"):
            rows = table.find_all("tr")
            if len(rows) < 2:
                continue

            # Find the header row
            header_idx, col_score, col_label, col_mutation = self._find_header(rows)
            if header_idx is None:
                continue

            # Locate the "No Mutation" row (original peptide)
            for row in rows[header_idx + 1 :]:
                cells = [td.get_text(strip=True) for td in row.find_all(["td", "th"])]
                if not cells:
                    continue

                # Accept the first data row OR the one that says "No Mutation"
                is_no_mutation = (
                    col_mutation is not None
                    and col_mutation < len(cells)
                    and "no mutation" in cells[col_mutation].lower()
                )
                is_first_data = row == rows[header_idx + 1]

                if not (is_no_mutation or is_first_data):
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
    ) -> tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
        """
        Search rows for a header containing score/prediction columns.

        Returns (header_row_idx, col_score, col_label, col_mutation).
        """
        for i, row in enumerate(rows):
            cells = [td.get_text(strip=True).lower() for td in row.find_all(["td", "th"])]

            # Must contain at least a score or prediction-like header
            has_score = any("score" in c for c in cells)
            has_pred = any("predict" in c or "result" in c for c in cells)
            if not (has_score or has_pred):
                continue

            col_score = next(
                (j for j, c in enumerate(cells) if "score" in c), None
            )
            col_label = next(
                (j for j, c in enumerate(cells) if "predict" in c or "result" in c),
                None,
            )
            col_mutation = next(
                (j for j, c in enumerate(cells) if "mutation" in c or "position" in c),
                None,
            )
            return i, col_score, col_label, col_mutation

        return None, None, None, None

    # ------------------------------------------------------------------
    # Regex helpers
    # ------------------------------------------------------------------

    def _regex_score(self, text: str) -> Optional[float]:
        for pattern in (self._SCORE_RE, self._SCORE_RE2):
            m = pattern.search(text)
            if m:
                return self._safe_float(m.group(1))
        return None

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
        """Standardize label variants from the server."""
        lower = label.lower()
        if "non" in lower:
            return "Non Anti-Angiogenic"
        return "Anti-Angiogenic"
