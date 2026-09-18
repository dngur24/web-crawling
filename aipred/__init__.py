"""
AIPpred Web Scraper Package
Anti-Inflammatory Peptide Prediction via http://www.thegleelab.org/AIPpred/
"""

from .predictor import Predictor
from .parser import Parser
from .models import PredictionResult

__all__ = ["Predictor", "Parser", "PredictionResult"]
__version__ = "1.0.0"
