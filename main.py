#!/usr/bin/env python3
"""
AntiAngioPred CLI / usage example.

Usage examples
--------------
# Single peptide (FullSeq method)
python main.py YCNINEVCHYARRNDKSYWL

# Single peptide with NT15 method
python main.py YCNINEVCHYARRNDKSYWL --method NT15

# Batch from a text file (one peptide per line)
python main.py --file peptides.txt --output results.csv

# Batch from command line
python main.py YCNINEVCHYARRNDKSYWL ACDEFGHIKLMNPQ --output results.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from antiangiopred import Predictor, Method, PredictionResult


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        level=level,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict anti-angiogenic activity of peptides via AntiAngioPred.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "peptides",
        nargs="*",
        metavar="SEQUENCE",
        help="One or more peptide sequences (amino-acid letters, 5–50 residues).",
    )
    parser.add_argument(
        "--file", "-f",
        type=Path,
        help="Path to a text file with one peptide per line.",
    )
    parser.add_argument(
        "--method", "-m",
        choices=["NT15", "FullSeq"],
        default="FullSeq",
        help="Prediction method (default: FullSeq).",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Path to save CSV output (optional).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Request timeout in seconds (default: 30).",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of retries on transient errors (default: 3).",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=1.0,
        help="Delay between batch requests in seconds (default: 1.0).",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose/debug logging.",
    )
    return parser.parse_args()


def collect_peptides(args: argparse.Namespace) -> list[str]:
    """Merge peptides from positional args and/or --file."""
    peptides: list[str] = list(args.peptides)

    if args.file:
        file_path: Path = args.file
        if not file_path.exists():
            print(f"[ERROR] File not found: {file_path}", file=sys.stderr)
            sys.exit(1)
        lines = file_path.read_text(encoding="utf-8").splitlines()
        peptides.extend(line.strip() for line in lines if line.strip())

    if not peptides:
        print("[ERROR] No peptide sequences provided. Use positional args or --file.")
        print("        Run with --help for usage details.")
        sys.exit(1)

    return peptides


def print_results(results: list[PredictionResult]) -> None:
    """Pretty-print results to stdout."""
    width = 52
    print("\n" + "=" * width)
    print(f"{'PEPTIDE':<22} {'METHOD':<8} {'LABEL':<25} {'SCORE'}")
    print("-" * width)
    for r in results:
        if r.success:
            score_str = f"{r.score:.4f}" if r.score is not None else "N/A"
            label_str = r.label or "N/A"
            print(f"{r.peptide:<22} {r.method:<8} {label_str:<25} {score_str}")
        else:
            print(f"{r.peptide:<22} {r.method:<8} ERROR: {r.error}")
    print("=" * width + "\n")


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)

    peptides = collect_peptides(args)

    predictor = Predictor(
        method=args.method,
        timeout=args.timeout,
        max_retries=args.retries,
        inter_request_delay=args.delay,
    )

    print(f"[AntiAngioPred] Method: {args.method} | Peptides: {len(peptides)}")

    if len(peptides) == 1:
        results = [predictor.predict(peptides[0])]
    else:
        results = predictor.predict_batch(peptides)

    print_results(results)

    if args.output:
        saved_path = PredictionResult.to_csv(results, args.output)
        print(f"[AntiAngioPred] Results saved to: {saved_path}")


if __name__ == "__main__":
    main()
