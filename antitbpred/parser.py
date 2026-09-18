"""
Parser: extracts structured data from the JS-rendered HTML returned by AntiTbPred.

Result page (disp.php?ran=XXXXX) — after Bootstrap-table JS execution — contains:
    ID | Seq | Score | Prediction | Steric hindrance | Amphipathicity | Net Hydrogen | pI

Key columns:
  - Score      : SVM raw score (float)
  - Prediction : classification result (float, threshold-based)

This parser handles:
 1. Table-based extraction (primary) — finds Score/Prediction columns
 2. Regex fallback on plain text
"""

from __future__ import annotations

import re
from typing import Optional

from bs4 import BeautifulSoup


class Parser:
    """
    Parses the JS-rendered HTML result page from AntiTbPred.

    The page contains a Bootstrap table with columns:
        ID | Seq | Score | Prediction | Steric hindrance | Amphipathicity | ...

    Returns (score, label) for the first valid data row.
    """

    # ── Regex fallbacks ──────────────────────────────────────────────────
    _SCORE_RE = re.compile(r"score\s*[:\-]?\s*([-\d.]+)", re.IGNORECASE)
    _LABEL_RE = re.compile(
        r"\b(AntiTB|Anti[- ]TB|Non[- ]AntiTB|Non[- ]Anti[- ]TB)\b",
        re.IGNORECASE,
    )

    def parse(self, html: str) -> tuple[Optional[float], Optional[str]]:
        """
        Extract (score, label) from the rendered result HTML.

        Parameters
        ----------
        html : str  – JS-rendered HTML of the result page

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
        return self._regex_score(text), self._regex_label(text)

    def parse_all(
        self, html: str
    ) -> list[tuple[str, Optional[float], Optional[str]]]:
        """
        Extract all result rows: list of (seq_id, score, label).
        """
        soup = BeautifulSoup(html, "lxml")
        results = []

        for table in soup.find_all("table"):
            rows = table.find_all("tr")
            if len(rows) < 2:
                continue

            header_idx, col_id, col_score, col_label = self._find_header(rows)
            if header_idx is None:
                continue

            for row in rows[header_idx + 1:]:
                cells = [td.get_text(strip=True) for td in row.find_all(["td", "th"])]
                if not cells or not any(c for c in cells if c and c != "-"):
                    continue

                seq_id = cells[col_id] if col_id is not None and col_id < len(cells) else ""
                score = (
                    self._safe_float(cells[col_score])
                    if col_score is not None and col_score < len(cells)
                    else None
                )
                label = self._derive_label(score, cells, col_label)
                results.append((seq_id, score, label))

        return results

    # ------------------------------------------------------------------
    # Table parsing
    # ------------------------------------------------------------------

    def _parse_table(
        self, soup: BeautifulSoup
    ) -> tuple[Optional[float], Optional[str]]:
        """Extract (score, label) from the first valid data row."""
        for table in soup.find_all("table"):
            rows = table.find_all("tr")
            if len(rows) < 2:
                continue

            header_idx, col_id, col_score, col_label = self._find_header(rows)
            if header_idx is None:
                continue

            for row in rows[header_idx + 1:]:
                cells = [td.get_text(strip=True) for td in row.find_all(["td", "th"])]
                if not cells or not any(c for c in cells if c and c != "-"):
                    continue

                score = (
                    self._safe_float(cells[col_score])
                    if col_score is not None and col_score < len(cells)
                    else None
                )
                label = self._derive_label(score, cells, col_label)
                if score is not None or label is not None:
                    return score, label

        return None, None

    def _find_header(
        self, rows: list
    ) -> tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
        """
        Locate header row with Score and Prediction columns.

        AntiTbPred columns: ID | Seq | Score | Prediction | ...

        Returns (header_row_idx, col_id, col_score, col_label).
        """
        for i, row in enumerate(rows):
            cells = [
                td.get_text(strip=True).lower()
                for td in row.find_all(["td", "th"])
            ]

            has_score = any("score" in c for c in cells)
            has_pred  = any("predict" in c or "label" in c for c in cells)
            if not (has_score or has_pred):
                continue

            col_id = next(
                (j for j, c in enumerate(cells) if c in ("id", "name", "seq_id")), None
            )
            col_score = next(
                (j for j, c in enumerate(cells) if "score" in c and "steric" not in c),
                None,
            )
            col_label = next(
                (j for j, c in enumerate(cells)
                 if "predict" in c or "label" in c or "class" in c),
                None,
            )
            return i, col_id, col_score, col_label

        return None, None, None, None

    # ------------------------------------------------------------------
    # Label derivation
    # ------------------------------------------------------------------

    def _derive_label(
        self,
        score: Optional[float],
        cells: list[str],
        col_label: Optional[int],
    ) -> Optional[str]:
        """
        Derive classification label for AntiTbPred results.

        The server's 'Prediction' column contains a numeric probability
        (not a text label). The SVM 'Score' column sign determines the class:
          - Score >= 0 → AntiTB
          - Score <  0 → Non-AntiTB

        If the Prediction column contains a recognizable text label
        (e.g. 'AntiTB', 'Non-AntiTB'), that takes priority.
        """
        # 1. Check if Prediction column contains a recognizable text label
        if col_label is not None and col_label < len(cells):
            cell_text = cells[col_label].strip()
            if cell_text and cell_text not in ("-", ""):
                # Try to parse as float — if it is, skip (it's a probability)
                try:
                    float(cell_text)
                except ValueError:
                    # It's text — normalize it
                    return self._normalize_label(cell_text)

        # 2. Derive from score sign (SVM convention: positive = predicted positive class)
        if score is not None:
            return "AntiTB" if score >= 0 else "Non-AntiTB"

        return None

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
            return float(value.strip())
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _normalize_label(label: str) -> str:
        """Standardize label variants from the server."""
        lower = label.strip().lower()
        if "non" in lower:
            return "Non-AntiTB"
        if "anti" in lower or "tb" in lower:
            return "AntiTB"
        return label.strip()
