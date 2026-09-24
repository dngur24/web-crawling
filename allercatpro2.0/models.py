"""
Data models for AllerCatPro 2.0 results.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Allergenicity evidence levels (AllerCatPro 2.0 terminology)
# ---------------------------------------------------------------------------

ALLERGENIC_LABELS: frozenset[str] = frozenset(
    {"strong evidence", "weak evidence"}
)
SAFE_LABEL: str = "no evidence"


@dataclass
class PredictionResult:
    """
    Standardized result from one AllerCatPro 2.0 prediction.

    Attributes
    ----------
    peptide : str
        Input peptide/protein sequence (upper-cased, stripped).
    name : str
        FASTA description line used as the query identifier.
    result : Optional[str]
        Evidence label returned by the server:
        'strong evidence', 'weak evidence', or 'no evidence'.
    comment : Optional[str]
        Short comment from the server (e.g. '6-mer x3 nolcr', 'no hits').
    sequence_length : Optional[int]
        Sequence length as reported by the server.
    gluten_repeats : Optional[int]
        Number of gluten-like Q-repeats detected.
    hexamer_overlaps : Optional[int]
        Number of 3×6-mer overlaps with known allergens.
    raw_html : Optional[str]
        Full HTML of the result page (not exported to CSV; useful for debug).
    error : Optional[str]
        Error message when the prediction failed.
    """

    peptide: str
    name: str = ""
    result: Optional[str] = None
    comment: Optional[str] = None
    sequence_length: Optional[int] = None
    gluten_repeats: Optional[int] = None
    hexamer_overlaps: Optional[int] = None
    raw_html: Optional[str] = field(default=None, repr=False)
    error: Optional[str] = None

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------

    @property
    def success(self) -> bool:
        """True when the prediction completed without errors."""
        return self.error is None

    @property
    def is_allergenic(self) -> bool:
        """True when the server reports any level of allergenic evidence."""
        if self.result is None:
            return False
        return self.result.lower() in ALLERGENIC_LABELS

    @property
    def is_safe(self) -> bool:
        """True when the server reports no allergenic evidence."""
        if self.result is None:
            return False
        return self.result.lower() == SAFE_LABEL

    def to_dict(self, include_raw_html: bool = False) -> dict:
        """Return a plain dictionary representation."""
        d = {
            "name": self.name,
            "peptide": self.peptide,
            "result": self.result,
            "comment": self.comment,
            "sequence_length": self.sequence_length,
            "gluten_repeats": self.gluten_repeats,
            "hexamer_overlaps": self.hexamer_overlaps,
            "is_allergenic": self.is_allergenic,
            "error": self.error,
        }
        if include_raw_html:
            d["raw_html"] = self.raw_html
        return d

    def to_json(self, indent: int = 2) -> str:
        """Serialize to a JSON string (without raw_html)."""
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=indent)

    # ------------------------------------------------------------------
    # Batch helpers (class-level)
    # ------------------------------------------------------------------

    @classmethod
    def to_csv(
        cls,
        results: list["PredictionResult"],
        path: str | Path,
        *,
        encoding: str = "utf-8-sig",  # BOM for Excel compatibility
    ) -> Path:
        """
        Write a list of PredictionResult objects to a CSV file.

        Parameters
        ----------
        results : list[PredictionResult]
        path : str or Path  – destination file path
        encoding : str

        Returns
        -------
        Path  – resolved path to the written file
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        fieldnames = [
            "name",
            "peptide",
            "result",
            "comment",
            "sequence_length",
            "gluten_repeats",
            "hexamer_overlaps",
            "is_allergenic",
            "error",
        ]
        with path.open("w", newline="", encoding=encoding) as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(r.to_dict())

        return path.resolve()
