"""
AntiCP 2.0 Web Scraper Package

Prediction of Anti-Cancer Peptides
Target: https://webs.iiitd.edu.in/raghava/anticp2/predict.php
"""

from .models import PredictionMethod, PredictionResult, ThresholdValue
from .parser import Parser
from .predictor import Predictor

__all__ = [
    "Predictor",
    "Parser",
    "PredictionResult",
    "PredictionMethod",
    "ThresholdValue",
]
__version__ = "1.0.0"
