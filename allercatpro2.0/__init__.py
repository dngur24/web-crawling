"""
AllerCatPro 2.0 Web Scraper Package

Predicts protein allergenicity potential via
https://allercatpro.bii.a-star.edu.sg/

Quick start
-----------
    from allercatpro2_0 import Predictor, PredictionResult

    predictor = Predictor()

    # Filter safe peptides (remove allergens)
    peptides = ["ACDEFGHIKLM", "MKTFLILALLA", ...]
    safe = predictor.filter_safe(peptides)

    # Full batch with CSV export
    results = predictor.predict_batch_bulk(peptides)
    PredictionResult.to_csv(results, "results.csv")
"""

from .predictor import Predictor
from .parser import Parser
from .models import PredictionResult, ALLERGENIC_LABELS, SAFE_LABEL

__all__ = [
    "Predictor",
    "Parser",
    "PredictionResult",
    "ALLERGENIC_LABELS",
    "SAFE_LABEL",
]
__version__ = "1.0.0"
