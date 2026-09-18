"""
AntiTbPred Web Scraper Package
Anti-Tubercular Peptide Prediction via https://webs.iiitd.edu.in/raghava/antitbpred/
"""

from .predictor import Predictor
from .parser import Parser
from .models import Method, PredictionResult

__all__ = ["Predictor", "Method", "Parser", "PredictionResult"]
__version__ = "1.0.0"
