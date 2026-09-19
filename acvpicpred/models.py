"""
Data models for ACVPICPred results.

Target: https://i.uestc.edu.cn/acvpICPred/main/Main.php
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional


# ──────────────────────────────────────────────────────────────────────────────
# Parameter Enums  (mirror <select> options from the HTML form)
# ──────────────────────────────────────────────────────────────────────────────

class VirusType(str, Enum):
    """Virus type options available on the ACVPICPred form."""
    FCoV        = "FCoV"
    HCoV_229E   = "HCoV-229E"
    HCoV_NL63   = "HCoV-NL63"
    HCoV_OC43   = "HCoV-OC43"
    MERS_CoV    = "MERS-CoV"
    SARS_CoV    = "SARS-CoV"
    SARS_CoV_2  = "SARS-CoV-2"
    TGEV        = "TGEV"
    OTHERS      = "others"


class AssayType(str, Enum):
    """Experiment / assay type options."""
    CELL_CELL_FUSION            = "cell-cell fusion"
    PLAQUE_REDUCTION            = "plaque reduction assay"
    CYTOPATHIC_EFFECT           = "cytopathic effect assay"
    PSEUDOTYPED_VIRUS           = "pseudotyped virus infection inhibition assay"
    BETA_GAL_FUSION             = "beta-gal complementation-based fusion assay"
    BIOTINYLATED_ELISA          = "biotinylated enzyme-linked immunosorbent assay"


class ValueType(str, Enum):
    """Inhibition / suppression value type."""
    IC50 = "IC50"
    EC50 = "EC50"
    IC90 = "IC90"


class UnitType(str, Enum):
    """Unit type for the inhibition value."""
    uM     = "uM"
    ug_mL  = "ug/ML"


# ──────────────────────────────────────────────────────────────────────────────
# Result dataclass
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class PredictionResult:
    """
    Standardized result from one ACVPICPred prediction.

    Attributes
    ----------
    peptide : str
        Input peptide sequence (upper-cased).
    pdb_path : str
        Path to the PDB structure file that was submitted.
    virus : str
        Virus type selected for the prediction.
    assay : str
        Experiment type selected.
    value_type : str
        Inhibition value type (IC50 / EC50 / IC90).
    unit : str
        Unit type (uM / ug/ML).
    predicted_value : Optional[float]
        Numerical predicted activity value returned by the server.
    label : Optional[str]
        Classification label returned by the server (e.g. "Active", "Inactive").
    raw_html : Optional[str]
        Full HTML of the result page (stored for debugging; not exported to CSV).
    error : Optional[str]
        Error message if the prediction failed, otherwise None.
    """

    peptide: str
    pdb_path: str
    virus: str
    assay: str
    value_type: str
    unit: str
    predicted_value: Optional[float] = None
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
            "peptide":         self.peptide,
            "pdb_path":        self.pdb_path,
            "virus":           self.virus,
            "assay":           self.assay,
            "value_type":      self.value_type,
            "unit":            self.unit,
            "predicted_value": self.predicted_value,
            "label":           self.label,
            "error":           self.error,
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
            "peptide", "pdb_path", "virus", "assay",
            "value_type", "unit", "predicted_value", "label", "error",
        ]
        with path.open("w", newline="", encoding=encoding) as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in results:
                writer.writerow(r.to_dict())

        return path.resolve()
