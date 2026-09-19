"""
ACVPICPred Web Scraper Package

Anti-Coronavirus Peptide Activity Predictor
Target: https://i.uestc.edu.cn/acvpICPred/main/Main.php
"""

from .models import PredictionResult, VirusType, AssayType, ValueType, UnitType
from .parser import Parser
from .predictor import Predictor

__all__ = [
    "Predictor",
    "Parser",
    "PredictionResult",
    "VirusType",
    "AssayType",
    "ValueType",
    "UnitType",
]
__version__ = "1.0.0"
