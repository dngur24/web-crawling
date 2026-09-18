"""
Data models for AagingPEPred results.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import csv
import json
from pathlib import Path


class Model(str):
    """
    Prediction model options available on the AagingPEPred server.

    Attributes
    ----------
    STOCHASTIC : str
        Stochastically Validated Model – recommended.
        Uses random feature selection across AA composition, dipeptides,
        and physicochemical properties.
    LITERATURE : str
        Literature-Validated Model.
        Trained on experimentally characterized proteins from SwissProt.
    TARGETED : str
        Targeted Feature Analysis Model.
        Focuses on statistically significant molecular descriptors.
    HIGH_PRECISION : str
        High-Precision Database Model.
        SwissProt-trained with computationally selected optimal features.
    """

    STOCHASTIC = "all_features_randomfeatures"
    LITERATURE = "all_features_swissprot"
    TARGETED = "rfe_features_randomfeatures"
    HIGH_PRECISION = "rfe_features_swissprot"

    # Human-readable aliases
    _ALIAS_MAP: dict[str, str] = {
        "stochastic": "all_features_randomfeatures",
        "all_features_randomfeatures": "all_features_randomfeatures",
        "literature": "all_features_swissprot",
        "all_features_swissprot": "all_features_swissprot",
        "targeted": "rfe_features_randomfeatures",
        "rfe_features_randomfeatures": "rfe_features_randomfeatures",
        "high_precision": "rfe_features_swissprot",
        "rfe_features_swissprot": "rfe_features_swissprot",
    }

    @classmethod
    def resolve(cls, value: str) -> str:
        """Resolve alias or raw value to a valid model string."""
        resolved = cls._ALIAS_MAP.get(value.lower())
        if resolved is None:
            raise ValueError(
                f"Unknown model '{value}'. "
                f"Valid values: {list(cls._ALIAS_MAP.keys())}"
            )
        return resolved


@dataclass
class PredictionResult:
    """
    Standardized result from one AagingPEPred prediction.

    Attributes
    ----------
    peptide : str
        Input peptide sequence (upper-cased).
    model : str
        Model used for prediction.
    score : Optional[float]
        Prediction score returned by the server (None if not parsed).
    label : Optional[str]
        Classification label, e.g. 'Anti-Aging' or 'Non Anti-Aging'.
    raw_html : Optional[str]
        Full HTML of the result page (for debugging; not exported to CSV).
    error : Optional[str]
        Error message if the prediction failed.
    """

    peptide: str
    model: str
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
            "model": self.model,
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

        fieldnames = ["peptide", "model", "score", "label", "error"]
        with path.open("w", newline="", encoding=encoding) as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(r.to_dict())

        return path.resolve()
