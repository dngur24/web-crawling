"""
Parser: extracts structured prediction data from the raw HTML returned
by AllerCatPro 2.0.

Result page table structure (from live observation):
    Row 0  – merged group headers (Protein sequence of interest, etc.)
    Row 1  – column headers: Query#, Name, Sequence Length, ..., Result, Comment
    Row 2+ – one data row per submitted sequence

    Header row has 18 cells, but data rows have 21 cells because three
    ``rowspan`` header columns ('Potential cross-reactivity',
    'Similarity to autoimmune allergen', 'Similarity to low allergenic
    protein') each contribute an extra ``<td>`` only in data rows.

    Fixed indices in DATA rows:
        [1]  Name
        [2]  Sequence Length
        [3]  Gluten allergens (# of Q-repeats)
        [4]  # of 3x6-mer overlaps
        [-4] Result  : 'strong evidence' | 'weak evidence' | 'no evidence'
        [-3] Comment : e.g. '6-mer x3 nolcr', 'no hits'
"""

from __future__ import annotations

import logging
from typing import Optional

from bs4 import BeautifulSoup

from .models import PredictionResult

logger = logging.getLogger(__name__)

# Fixed column indices in DATA rows (not header rows)
_IDX_NAME    = 1
_IDX_SEQ_LEN = 2
_IDX_GLUTEN  = 3
_IDX_HEXAMER = 4
_IDX_RESULT  = -4   # 4th from the end
_IDX_COMMENT = -3   # 3rd from the end


class Parser:
    """
    Parses the raw HTML result page from AllerCatPro 2.0.

    The server returns a single-page result containing a ``<table
    class="tableSMS">`` with two header rows and one data row per query
    sequence.
    """

    def parse(self, html: str) -> list[PredictionResult]:
        """
        Extract a list of :class:`PredictionResult` from the result HTML.

        Parameters
        ----------
        html : str
            Raw HTML of the AllerCatPro 2.0 result page.

        Returns
        -------
        list[PredictionResult]
            One entry per submitted sequence, in submission order.
            Returns an empty list when parsing fails entirely.
        """
        soup = BeautifulSoup(html, "lxml")

        # Error page detection
        error_h1 = soup.find("h1")
        if error_h1:
            h1_text = error_h1.get_text().lower()
            if "not in fasta" in h1_text or "too large" in h1_text:
                logger.error(
                    "Server returned an error page: %s",
                    error_h1.get_text(strip=True),
                )
                return []

        main_table = soup.find("table", class_="tableSMS")
        if main_table is None:
            logger.warning("Could not find <table class='tableSMS'> in response.")
            return []

        return self._parse_table(main_table)

    # ------------------------------------------------------------------
    # Table parsing
    # ------------------------------------------------------------------

    def _parse_table(self, table) -> list[PredictionResult]:
        rows = table.find_all("tr")
        if len(rows) < 3:
            logger.warning("Result table has fewer than 3 rows; cannot parse.")
            return []

        # Row 0: merged group headers  (skip)
        # Row 1: per-column headers    (skip – we use fixed indices)
        # Row 2+: data rows
        results: list[PredictionResult] = []
        for row in rows[2:]:
            cells = [td.get_text(strip=True) for td in row.find_all(["td", "th"])]
            if not cells:
                continue
            result = self._parse_data_row(cells)
            if result is not None:
                results.append(result)

        return results

    @staticmethod
    def _parse_data_row(cells: list[str]) -> Optional[PredictionResult]:
        """Parse a single data row into a :class:`PredictionResult`."""
        def get(idx: int) -> Optional[str]:
            try:
                v = cells[idx].strip()
                return v if v and v != "-" else None
            except IndexError:
                return None

        name       = get(_IDX_NAME) or ""
        result_raw = get(_IDX_RESULT)
        comment    = get(_IDX_COMMENT)
        seq_len    = _safe_int(get(_IDX_SEQ_LEN))
        gluten     = _safe_int(get(_IDX_GLUTEN))
        hexamer    = _safe_int(get(_IDX_HEXAMER))

        if result_raw is None:
            logger.warning(
                "No 'result' value at index %d in row (len=%d): %s",
                _IDX_RESULT,
                len(cells),
                cells,
            )
            return None

        return PredictionResult(
            peptide="",          # filled by Predictor after submission
            name=name,
            result=result_raw.lower(),
            comment=comment,
            sequence_length=seq_len,
            gluten_repeats=gluten,
            hexamer_overlaps=hexamer,
        )


# ---------------------------------------------------------------------------
# Module-level helper
# ---------------------------------------------------------------------------

def _safe_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None
