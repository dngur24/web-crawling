"""
AagingPEPred Web Scraper Package
Anti-Aging Peptide Prediction via https://project.iith.ac.in/cgntlab/aagingpepred/
"""

from .predictor import Predictor
from .parser import Parser
from .models import Model, PredictionResult

__all__ = ["Predictor", "Model", "Parser", "PredictionResult"]
__version__ = "1.0.0"
