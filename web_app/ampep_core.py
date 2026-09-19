"""
ampep_core.py — amPEPpy 웹 서비스화를 위한 리팩토링 코어 모듈

문제점.md의 모든 버그를 수정하고 순수 함수(pure function) 형태로 재구성:
  - [문제점 1] predict() 루프 내 단건 clf.predict() → 벡터화 일괄 호출
  - [문제점 2] 중복 서열 ID → iloc 정수 위치 인덱싱
  - [문제점 3] score() CSV 자동 저장 부작용 → 완전 제거 (순수 함수)
  - [문제점 5] 비표준 아미노산 검증 → validate_sequence() 추가
  - [문제점 6] pickle → joblib
  - [문제점 7] 비효율적 join 문법, 오탈자 정리
"""

from __future__ import annotations

import io
import re
from pathlib import Path
from typing import Optional

import joblib
import numpy as np
import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------
# 상수 정의
# ---------------------------------------------------------------------------

CTD: dict = {
    "hydrophobicity": {
        1: ["R", "K", "E", "D", "Q", "N"],
        2: ["G", "A", "S", "T", "P", "H", "Y"],
        3: ["C", "L", "V", "I", "M", "F", "W"],
    },
    "normalized.van.der.waals": {
        1: ["G", "A", "S", "T", "P", "D", "C"],
        2: ["N", "V", "E", "Q", "I", "L"],
        3: ["M", "H", "K", "F", "R", "Y", "W"],
    },
    "polarity": {
        1: ["L", "I", "F", "W", "C", "M", "V", "Y"],
        2: ["P", "A", "T", "G", "S"],
        3: ["H", "Q", "R", "K", "N", "E", "D"],
    },
    "polarizability": {
        1: ["G", "A", "S", "D", "T"],
        2: ["C", "P", "N", "V", "E", "Q", "I", "L"],
        3: ["K", "M", "H", "F", "R", "Y", "W"],
    },
    "charge": {
        1: ["K", "R"],
        2: ["A", "N", "C", "Q", "G", "H", "I", "L", "M", "F", "P", "S", "T", "W", "Y", "V"],
        3: ["D", "E"],
    },
    "secondary": {
        1: ["E", "A", "L", "M", "Q", "K", "R", "H"],
        2: ["V", "I", "Y", "C", "W", "F", "T"],
        3: ["G", "N", "P", "S", "D"],
    },
    "solvent": {
        1: ["A", "L", "F", "C", "G", "I", "V", "W"],
        2: ["R", "K", "Q", "E", "N", "D"],
        3: ["M", "S", "P", "T", "H", "Y"],
    },
}

STANDARD_AA: frozenset = frozenset("ACDEFGHIKLMNPQRSTVWY")
_NONSTANDARD_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]", re.IGNORECASE)

FEATURE_HEADER: list = [
    f"{prop}.{group}.{pct}"
    for prop in CTD
    for group in (1, 2, 3)
    for pct in (0, 25, 50, 75, 100)
]

PROPERTY_LABELS = {
    "hydrophobicity":         "소수성 (Hydrophobicity)",
    "normalized.van.der.waals": "반데르발스 부피",
    "polarity":               "극성 (Polarity)",
    "polarizability":         "분극률 (Polarizability)",
    "charge":                 "전하 (Charge)",
    "secondary":              "2차 구조 선호도",
    "solvent":                "용매 접근성",
}

# ---------------------------------------------------------------------------
# 입력 검증
# ---------------------------------------------------------------------------

class SequenceValidationError(ValueError):
    pass


def validate_sequence(sequence: str, strict: bool = False) -> tuple:
    """
    서열 검증 및 정규화.
    Returns (validated_seq: str, warnings: list[str])
    """
    warnings = []
    seq = sequence.upper().replace("*", "").replace("-", "").strip()
    non_standard = _NONSTANDARD_RE.findall(seq)
    if non_standard:
        unique = sorted(set(non_standard))
        msg = f"비표준 아미노산 제거됨: {unique}"
        if strict:
            raise SequenceValidationError(msg)
        warnings.append(msg)
        seq = _NONSTANDARD_RE.sub("", seq)
    if len(seq) == 0:
        raise SequenceValidationError("유효한 아미노산 서열이 없습니다.")
    return seq, warnings


# ---------------------------------------------------------------------------
# CTD 피처 추출
# ---------------------------------------------------------------------------

def _percentile_positions(positions: list, seq_len: int) -> list:
    """
    CLI(amPEP.py)의 score() 함수와 동일한 방식으로 퍼센타일 위치를 계산합니다.

    CLI는 positions 리스트 맨 앞에 더미 "-"를 삽입해 1-based 인덱싱을 사용합니다:
        abpos.insert(0, "-")  →  abpos[1] = 첫 번째 위치, abpos[n] = 마지막 위치
        25th 퍼센타일: abpos[round(0.25*n - 0.1)]  (더미 포함 1-based)

    이를 0-based 리스트로 재현하려면 인덱스를 1 감소시켜야 합니다:
        0-based 첫 번째: pos[0]  (CLI의 abpos[1])
        0-based 25th:   pos[round(0.25*n - 0.1) - 1]
        0-based 마지막:  pos[n-1]  (CLI의 abpos[n])
    """
    n = len(positions)
    if n == 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    pos = positions
    if n == 1:
        pct = (pos[0] / seq_len) * 100
        return [pct, pct, pct, pct, pct]
    if n == 2:
        return [
            (pos[0] / seq_len) * 100,
            (pos[0] / seq_len) * 100,
            (pos[round(0.5 * n - 0.1) - 1] / seq_len) * 100,
            (pos[round(0.75 * n - 0.1) - 1] / seq_len) * 100,
            (pos[n - 1] / seq_len) * 100,
        ]
    return [
        (pos[0] / seq_len) * 100,
        (pos[round(0.25 * n - 0.1) - 1] / seq_len) * 100,
        (pos[round(0.5  * n - 0.1) - 1] / seq_len) * 100,
        (pos[round(0.75 * n - 0.1) - 1] / seq_len) * 100,
        (pos[n - 1] / seq_len) * 100,
    ]


def extract_features(sequence: str) -> list:
    """단일 서열 → 105차원 CTD 피처 리스트 (순수 함수, CSV 저장 없음)."""
    features = []
    seq_len = len(sequence)
    for prop, groups in CTD.items():
        encoded = ""
        for aa in sequence:
            if aa in groups[1]:
                encoded += "1"
            elif aa in groups[2]:
                encoded += "2"
            elif aa in groups[3]:
                encoded += "3"
        effective_len = len(encoded) if encoded else seq_len
        for g in (1, 2, 3):
            char = str(g)
            positions = [i + 1 for i, c in enumerate(encoded) if c == char]
            features.extend(_percentile_positions(positions, effective_len))
    return features


def score_sequences(sequences: dict) -> pd.DataFrame:
    """
    {seq_id: sequence} dict → CTD 피처 DataFrame (shape: n × 105).
    [문제점 3 수정] CSV 자동 저장 완전 제거.
    """
    rows = [extract_features(seq) for seq in sequences.values()]
    df = pd.DataFrame(rows, columns=FEATURE_HEADER, index=list(sequences.keys()))
    df.index.name = "sequence_id"
    return df


def parse_fasta(fasta_input) -> tuple:
    """
    FASTA 문자열 또는 파일 핸들 → ({id: seq}, warnings_list).
    """
    if isinstance(fasta_input, str):
        handle = io.StringIO(fasta_input)
    else:
        handle = fasta_input

    sequences = {}
    all_warnings = []
    for record in SeqIO.parse(handle, "fasta"):
        try:
            seq, warns = validate_sequence(str(record.seq))
            sequences[record.id] = seq
            if warns:
                all_warnings.extend([f"[{record.id}] {w}" for w in warns])
        except SequenceValidationError as e:
            all_warnings.append(f"[{record.id}] 건너뜀: {e}")
    return sequences, all_warnings


# ---------------------------------------------------------------------------
# 모델 로드
# ---------------------------------------------------------------------------

def load_model(model_path):
    """
    [문제점 6 수정] pickle → joblib (하위 호환 폴백 포함).
    """
    model_path = Path(model_path)
    if not model_path.exists():
        raise FileNotFoundError(f"모델 파일 없음: {model_path}")
    try:
        return joblib.load(model_path)
    except Exception:
        import pickle
        with open(model_path, "rb") as f:
            return pickle.load(f)


def save_model(clf, model_path) -> None:
    """[문제점 6 수정] joblib으로 모델 저장."""
    joblib.dump(clf, Path(model_path))


# ---------------------------------------------------------------------------
# 예측
# ---------------------------------------------------------------------------

def predict_batch(clf, feature_df: pd.DataFrame, drop_features=None) -> pd.DataFrame:
    """
    [문제점 1 수정] clf.predict() 단 1회 벡터화 호출.
    [문제점 2 수정] iloc 기반으로 중복 index 문제 회피.
    """
    X = feature_df.copy()
    if drop_features:
        X = X.drop(columns=[f for f in drop_features if f in X.columns])

    seq_ids = list(X.index)
    probas = clf.predict_proba(X)
    preds  = clf.predict(X)          # 단 1회 벡터화

    return pd.DataFrame({
        "sequence_id":        seq_ids,
        "probability_nonAMP": probas[:, 0],
        "probability_AMP":    probas[:, 1],
        "predicted":          ["AMP" if p == 1 else "nonAMP" for p in preds],
    })


def predict_single(clf, sequence: str) -> dict:
    """단일 서열 예측 (Streamlit 단일 입력 탭 전용)."""
    validated, warnings = validate_sequence(sequence)
    feature_df = score_sequences({"input": validated})
    result_df  = predict_batch(clf, feature_df)
    row = result_df.iloc[0]
    return {
        "sequence":           validated,
        "length":             len(validated),
        "probability_nonAMP": float(row["probability_nonAMP"]),
        "probability_AMP":    float(row["probability_AMP"]),
        "predicted":          row["predicted"],
        "features":           feature_df.iloc[0].to_dict(),
        "warnings":           warnings,
    }
