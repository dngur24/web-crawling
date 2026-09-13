"""
Data models for AntiAngioPred results.
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Optional
import csv
import json
from pathlib import Path


@dataclass
class PredictionResult:
    """
    Standardized result from one AntiAngioPred prediction.

    Attributes
    ----------
    peptide : str
        Input peptide sequence (upper-cased).
    method : str
        Prediction method used: 'NT15' or 'FullSeq'.
    score : Optional[float]
        Raw prediction score returned by the server (None if not parsed).
    label : Optional[str]
        Classification label, e.g. 'Anti-Angiogenic' or 'Non Anti-Angiogenic'.
    raw_html : Optional[str]
        Full HTML of the result page (stored for debugging, not exported to CSV).
    error : Optional[str]
        Error message if the prediction failed.
    """

    peptide: str
    method: str
    score: Optional[float] = None
    label: Optional[str] = None
    raw_html: Optional[str] = field(default=None, repr=False)
    error: Optional[str] = None

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------

    @property
    def success(self) -> bool:
        """True when the prediction completed without errors."""
        return self.error is None

    def to_dict(self, include_raw_html: bool = False) -> dict:
        """Return a plain dictionary representation."""
        d = {
            "peptide": self.peptide,
            "method": self.method,
            "score": self.score,
            "label": self.label,
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

        fieldnames = ["peptide", "method", "score", "label", "error"]
        with path.open("w", newline="", encoding=encoding) as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(r.to_dict())

        return path.resolve()
