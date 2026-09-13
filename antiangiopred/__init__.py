"""
AntiAngioPred Web Scraper Package
Anti-Angiogenic Peptide Prediction via https://webs.iiitd.edu.in/raghava/antiangiopred/
"""

from .predictor import Predictor, Method
from .parser import Parser
from .models import PredictionResult

__all__ = ["Predictor", "Method", "Parser", "PredictionResult"]
__version__ = "1.0.0"
