"""
app.py — amPEPpy Streamlit 웹 애플리케이션

탭 구성:
  1. 단일 서열 분석   — 즉석 AMP 예측 + 물리화학적 속성 레이더 + 서열 하이라이팅
  2. 배치 분석 (FASTA) — 다중 서열 업로드, 결과 테이블, CTD 히트맵, TSV 다운로드
"""

import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# 현재 디렉터리를 경로에 추가 (ampep_core import용)
sys.path.insert(0, str(Path(__file__).parent))
import ampep_core as core

# ---------------------------------------------------------------------------
# 페이지 설정
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="amPEPpy — 인실리코 항균 펩타이드 분석",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------------------------------------------------------------------------
# 모델 로드 (서버 시작 시 1회, 캐시)
# ---------------------------------------------------------------------------
@st.cache_resource(show_spinner="🔄 모델 로딩 중...")
def get_model(model_path: str):
    return core.load_model(model_path)


# ---------------------------------------------------------------------------
# 사이드바
# ---------------------------------------------------------------------------
with st.sidebar:
    st.image("https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/Peptide_bond.svg/320px-Peptide_bond.svg.png",
             width=160)
    st.title("⚙️ 설정")

    # 모델 경로 고정: pretrained_models/amPEP.model
    model_path = str(Path(__file__).parent.parent / "pretrained_models" / "amPEP.model")
    if Path(model_path).exists():
        st.success(f"✅ 모델: `{Path(model_path).name}`")
    else:
        st.error(f"❌ 모델 파일을 찾을 수 없습니다: `{model_path}`")
        st.stop()

    st.markdown("---")
    st.markdown(
        "**amPEPpy**는 CTD Distribution 피처(105차원)와 "
        "Random Forest 분류기를 사용하여 항균 펩타이드(AMP)를 예측합니다.\n\n"
        "[GitHub](https://github.com/tlawrence3/amPEPpy)"
    )

# 모델 로드
clf = get_model(model_path)

# ---------------------------------------------------------------------------
# 헤더
# ---------------------------------------------------------------------------
st.title("🧬 amPEPpy: 인실리코 항균 펩타이드 분석")
st.caption("Random Forest + CTD Distribution 기반 AMP 예측 웹 서비스")

tab1, tab2 = st.tabs(["🔬 단일 서열 분석", "📂 배치 분석 (FASTA)"])


# ===========================================================================
# 탭 1: 단일 서열 분석
# ===========================================================================
with tab1:
    st.subheader("단일 아미노산 서열 예측")

    sequence_input = st.text_area(
        "아미노산 서열을 입력하세요 (1문자 코드)",
        height=100,
        placeholder="예: GIGKFLKKAKKFGKAFVKILKK",
        key="single_seq",
    )

    predict_btn = st.button("🚀 예측 실행", type="primary", key="predict_single")

    if predict_btn and sequence_input.strip():
        with st.spinner("분석 중..."):
            try:
                result = core.predict_single(clf, sequence_input.strip())
            except core.SequenceValidationError as e:
                st.error(f"❌ 서열 오류: {e}")
                st.stop()

        # 경고 표시
        for w in result["warnings"]:
            st.warning(f"⚠️ {w}")

        # ── 예측 결과 카드 ──────────────────────────────────────────────
        amp_prob    = result["probability_AMP"]
        nonamp_prob = result["probability_nonAMP"]
        predicted   = result["predicted"]
        seq_len     = result["length"]

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            label = "🟢 AMP" if predicted == "AMP" else "🔴 비-AMP"
            st.metric("예측 결과", label)
        with col2:
            st.metric("AMP 확률", f"{amp_prob:.3f}", delta=f"{amp_prob-0.5:+.3f} (기준 0.5)")
        with col3:
            st.metric("비-AMP 확률", f"{nonamp_prob:.3f}")
        with col4:
            st.metric("서열 길이", f"{seq_len} aa")

        # 확률 게이지
        fig_gauge = go.Figure(go.Indicator(
            mode="gauge+number",
            value=amp_prob * 100,
            domain={"x": [0, 1], "y": [0, 1]},
            title={"text": "AMP 확률 (%)"},
            gauge={
                "axis": {"range": [0, 100]},
                "bar":  {"color": "#2ecc71" if predicted == "AMP" else "#e74c3c"},
                "steps": [
                    {"range": [0,  50], "color": "#fadbd8"},
                    {"range": [50, 100], "color": "#d5f5e3"},
                ],
                "threshold": {
                    "line": {"color": "black", "width": 3},
                    "thickness": 0.75,
                    "value": 50,
                },
            },
        ))
        fig_gauge.update_layout(height=280, margin=dict(t=40, b=10, l=10, r=10))
        st.plotly_chart(fig_gauge, use_container_width=True)

        # ── 서열 잔기 하이라이팅 ──────────────────────────────────────
        st.markdown("#### 🎨 서열 물리화학적 속성 하이라이팅")
        prop_choice = st.selectbox(
            "표시할 속성",
            options=list(core.PROPERTY_LABELS.keys()),
            format_func=lambda k: core.PROPERTY_LABELS[k],
            key="prop_hl",
        )
        groups = core.CTD[prop_choice]
        colors = {1: "#3498db", 2: "#e67e22", 3: "#2ecc71"}
        group_labels = {1: "그룹 1 (저)", 2: "그룹 2 (중)", 3: "그룹 3 (고)"}

        html_seq = ""
        for aa in result["sequence"]:
            grp = next((g for g, aas in groups.items() if aa in aas), None)
            color = colors.get(grp, "#aaaaaa")
            html_seq += (
                f'<span style="background:{color};color:white;'
                f'padding:1px 4px;margin:1px;border-radius:3px;'
                f'font-family:monospace;font-size:1rem;">{aa}</span>'
            )
        st.markdown(html_seq, unsafe_allow_html=True)

        # 범례
        legend_html = " ".join(
            f'<span style="background:{colors[g]};color:white;padding:2px 8px;border-radius:3px;">'
            f'{group_labels[g]}</span>'
            for g in (1, 2, 3)
        )
        st.markdown(legend_html + ' <span style="color:#aaa">■ 비해당</span>', unsafe_allow_html=True)

        # ── 레이더 차트 (속성별 평균 CTD 점수) ──────────────────────
        st.markdown("#### 📡 물리화학적 속성 분포 레이더")
        features = result["features"]
        radar_vals = []
        radar_labels = []
        for prop, label in core.PROPERTY_LABELS.items():
            prop_feats = [v for k, v in features.items() if k.startswith(prop)]
            radar_vals.append(np.mean(prop_feats))
            radar_labels.append(label)

        fig_radar = go.Figure(go.Scatterpolar(
            r=radar_vals + [radar_vals[0]],
            theta=radar_labels + [radar_labels[0]],
            fill="toself",
            fillcolor="rgba(52,152,219,0.25)",
            line_color="#2980b9",
        ))
        fig_radar.update_layout(
            polar=dict(radialaxis=dict(range=[0, 100])),
            height=380,
            margin=dict(t=40, b=40, l=60, r=60),
        )
        st.plotly_chart(fig_radar, use_container_width=True)

        # ── CTD 피처 상세 히트맵 ────────────────────────────────────
        with st.expander("📊 CTD 피처 상세 보기 (105차원)"):
            feat_df = pd.DataFrame(
                [list(features.values())],
                columns=list(features.keys()),
            ).T.rename(columns={0: "값"})
            feat_df["속성"] = feat_df.index.str.rsplit(".", n=2).str[0]
            feat_df["그룹"] = feat_df.index.str.rsplit(".", n=2).str[1]
            feat_df["백분위"] = feat_df.index.str.rsplit(".", n=2).str[2]

            fig_hm = px.imshow(
                feat_df["값"].values.reshape(7, 15),
                x=[f"G{g}.{p}" for g in ["1","2","3"] for p in ["0","25","50","75","100"]],
                y=list(core.PROPERTY_LABELS.values()),
                color_continuous_scale="RdYlGn",
                zmin=0, zmax=100,
                aspect="auto",
                title="CTD Distribution 히트맵",
            )
            fig_hm.update_layout(height=320, margin=dict(t=40, b=10))
            st.plotly_chart(fig_hm, use_container_width=True)

    elif predict_btn:
        st.warning("⚠️ 서열을 입력해주세요.")


# ===========================================================================
# 탭 2: 배치 분석
# ===========================================================================
with tab2:
    st.subheader("FASTA 파일 배치 분석")

    col_up, col_paste = st.columns(2)
    with col_up:
        uploaded_fasta = st.file_uploader("FASTA 파일 업로드", type=["fasta", "fa", "faa", "txt"])
    with col_paste:
        pasted_fasta = st.text_area(
            "또는 FASTA 텍스트 직접 붙여넣기",
            height=150,
            placeholder=">seq1\nGIGKFLKKAKKFGKAFVKILKK\n>seq2\nMELITTIN...",
        )

    run_batch = st.button("🚀 배치 예측 실행", type="primary", key="predict_batch")

    # ── 예측 실행 (버튼 클릭 시에만) ────────────────────────────────
    if run_batch:
        if uploaded_fasta:
            fasta_text = uploaded_fasta.read().decode("utf-8", errors="replace")
        elif pasted_fasta.strip():
            fasta_text = pasted_fasta.strip()
        else:
            st.warning("⚠️ FASTA 파일을 업로드하거나 서열을 붙여넣으세요.")
            st.stop()

        with st.spinner("서열 파싱 및 피처 추출 중..."):
            sequences, parse_warnings = core.parse_fasta(fasta_text)

        if not sequences:
            st.error("❌ 유효한 서열이 없습니다. FASTA 형식을 확인하세요.")
            st.stop()

        with st.spinner("RF 모델 예측 중..."):
            feature_df = core.score_sequences(sequences)
            result_df  = core.predict_batch(clf, feature_df)

        # 결과를 session_state에 저장
        st.session_state["batch_result_df"]  = result_df
        st.session_state["batch_feature_df"] = feature_df
        st.session_state["batch_warnings"]   = parse_warnings

    # ── 저장된 결과가 있으면 표시 (정렬 변경 시에도 재사용) ─────────
    if "batch_result_df" in st.session_state:
        result_df   = st.session_state["batch_result_df"]
        feature_df  = st.session_state["batch_feature_df"]
        parse_warnings = st.session_state["batch_warnings"]

        # 경고 표시
        if parse_warnings:
            with st.expander(f"⚠️ 파싱 경고 {len(parse_warnings)}건"):
                for w in parse_warnings:
                    st.warning(w)

        st.info(f"✅ {len(result_df)}개 서열 예측 완료")

        # ── 요약 통계 ───────────────────────────────────────────────
        n_amp    = (result_df["predicted"] == "AMP").sum()
        n_nonamp = (result_df["predicted"] == "nonAMP").sum()
        mean_prob = result_df["probability_AMP"].mean()

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("총 서열 수",    len(result_df))
        c2.metric("AMP 예측",      n_amp,    delta=f"{n_amp/len(result_df)*100:.1f}%")
        c3.metric("비-AMP 예측",   n_nonamp, delta=f"{n_nonamp/len(result_df)*100:.1f}%")
        c4.metric("평균 AMP 확률", f"{mean_prob:.3f}")

        # ── AMP 확률 분포 히스토그램 ─────────────────────────────────
        fig_hist = px.histogram(
            result_df,
            x="probability_AMP",
            color="predicted",
            nbins=50,
            color_discrete_map={"AMP": "#2ecc71", "nonAMP": "#e74c3c"},
            title="AMP 확률 분포",
            labels={"probability_AMP": "AMP 확률", "count": "서열 수"},
        )
        fig_hist.add_vline(x=0.5, line_dash="dash", line_color="black", annotation_text="임계값 0.5")
        fig_hist.update_layout(height=320, margin=dict(t=40, b=10))
        st.plotly_chart(fig_hist, use_container_width=True)

        # ── 결과 테이블 ─────────────────────────────────────────────
        st.markdown("#### 📋 예측 결과 테이블")
        display_df = result_df.copy()
        display_df["probability_AMP"]    = display_df["probability_AMP"].round(4)
        display_df["probability_nonAMP"] = display_df["probability_nonAMP"].round(4)

        def highlight_amp(row):
            color = "#d5f5e3" if row["predicted"] == "AMP" else "#fadbd8"
            return [f"background-color: {color}"] * len(row)

        st.dataframe(
            display_df.style.apply(highlight_amp, axis=1),
            use_container_width=True,
            height=400,
        )

        # ── AMP / non-AMP 분리 표 ───────────────────────────────────
        st.markdown("#### 🔎 AMP / non-AMP 분리 결과")

        _sort_cols = ["sequence_id", "probability_AMP", "probability_nonAMP"]

        amp_base    = display_df[display_df["predicted"] == "AMP"].reset_index(drop=True)
        nonamp_base = display_df[display_df["predicted"] == "nonAMP"].reset_index(drop=True)

        split_tab_amp, split_tab_nonamp = st.tabs([
            f"🟢 AMP ({len(amp_base)}개)",
            f"🔴 non-AMP ({len(nonamp_base)}개)",
        ])

        with split_tab_amp:
            if amp_base.empty:
                st.info("AMP로 예측된 서열이 없습니다.")
            else:
                sort_col_amp, sort_dir_amp = st.columns([2, 1])
                with sort_col_amp:
                    amp_sort_by = st.selectbox(
                        "정렬 기준",
                        options=_sort_cols,
                        index=1,
                        key="amp_sort_col",
                    )
                with sort_dir_amp:
                    amp_asc = st.radio(
                        "정렬 방향",
                        options=["⬆ 오름차순", "⬇ 내림차순"],
                        index=1,
                        key="amp_sort_dir",
                        horizontal=True,
                    ) == "⬆ 오름차순"

                amp_df = amp_base.sort_values(amp_sort_by, ascending=amp_asc).reset_index(drop=True)

                def _hl_amp(row):
                    return ["background-color: #d5f5e3"] * len(row)
                st.dataframe(
                    amp_df.style.apply(_hl_amp, axis=1),
                    use_container_width=True,
                    height=min(400, max(120, len(amp_df) * 38 + 40)),
                )

        with split_tab_nonamp:
            if nonamp_base.empty:
                st.info("non-AMP로 예측된 서열이 없습니다.")
            else:
                sort_col_nonamp, sort_dir_nonamp = st.columns([2, 1])
                with sort_col_nonamp:
                    nonamp_sort_by = st.selectbox(
                        "정렬 기준",
                        options=_sort_cols,
                        index=1,
                        key="nonamp_sort_col",
                    )
                with sort_dir_nonamp:
                    nonamp_asc = st.radio(
                        "정렬 방향",
                        options=["⬆ 오름차순", "⬇ 내림차순"],
                        index=0,
                        key="nonamp_sort_dir",
                        horizontal=True,
                    ) == "⬆ 오름차순"

                nonamp_df = nonamp_base.sort_values(nonamp_sort_by, ascending=nonamp_asc).reset_index(drop=True)

                def _hl_nonamp(row):
                    return ["background-color: #fadbd8"] * len(row)
                st.dataframe(
                    nonamp_df.style.apply(_hl_nonamp, axis=1),
                    use_container_width=True,
                    height=min(400, max(120, len(nonamp_df) * 38 + 40)),
                )

        # ── CTD 히트맵 (상위 50개) ──────────────────────────────────
        if len(feature_df) <= 200:
            with st.expander("🗺️ CTD 피처 히트맵 (서열 간 비교)"):
                n_show = min(50, len(feature_df))
                hm_data = feature_df.iloc[:n_show]
                fig_batch_hm = px.imshow(
                    hm_data.values,
                    x=list(feature_df.columns),
                    y=list(hm_data.index),
                    color_continuous_scale="RdYlGn",
                    zmin=0, zmax=100,
                    aspect="auto",
                    title=f"CTD Distribution 히트맵 (상위 {n_show}개 서열)",
                )
                fig_batch_hm.update_layout(height=max(300, n_show * 14))
                st.plotly_chart(fig_batch_hm, use_container_width=True)

        # ── 다운로드 ────────────────────────────────────────────────
        st.markdown("#### 💾 결과 다운로드")

        # 필터링용 ID 집합
        amp_ids    = set(result_df.loc[result_df["predicted"] == "AMP",    "sequence_id"])
        nonamp_ids = set(result_df.loc[result_df["predicted"] == "nonAMP", "sequence_id"])

        amp_result_df    = result_df[result_df["sequence_id"].isin(amp_ids)].reset_index(drop=True)
        nonamp_result_df = result_df[result_df["sequence_id"].isin(nonamp_ids)].reset_index(drop=True)

        amp_feat_df    = feature_df[feature_df.index.isin(amp_ids)]
        nonamp_feat_df = feature_df[feature_df.index.isin(nonamp_ids)]

        def _to_tsv(df):
            buf = io.StringIO()
            df.to_csv(buf, sep="\t", index=False)
            return buf.getvalue()

        def _to_csv(df):
            buf = io.StringIO()
            df.to_csv(buf, index=True)
            return buf.getvalue()

        # 헤더 행
        _, hcol1, hcol2 = st.columns([1.2, 2, 2])
        hcol1.markdown("**예측 결과 (TSV)**")
        hcol2.markdown("**CTD 피처 (CSV)**")

        st.markdown("<hr style='margin:4px 0'>", unsafe_allow_html=True)

        # ── 전체 ──
        lbl_all, col_all_tsv, col_all_csv = st.columns([1.2, 2, 2])
        lbl_all.markdown("🗂️ **전체**")
        with col_all_tsv:
            st.download_button(
                "📥 전체 예측 결과.tsv",
                data=_to_tsv(result_df),
                file_name="ampep_all_predictions.tsv",
                mime="text/tab-separated-values",
                use_container_width=True,
                key="dl_all_tsv",
            )
        with col_all_csv:
            st.download_button(
                "📥 전체 CTD 피처.csv",
                data=_to_csv(feature_df),
                file_name="ampep_all_features.csv",
                mime="text/csv",
                use_container_width=True,
                key="dl_all_csv",
            )

        # ── AMP ──
        lbl_amp, col_amp_tsv, col_amp_csv = st.columns([1.2, 2, 2])
        lbl_amp.markdown("🟢 **AMP**")
        with col_amp_tsv:
            st.download_button(
                f"📥 AMP 예측 결과.tsv  ({len(amp_result_df)}개)",
                data=_to_tsv(amp_result_df),
                file_name="ampep_AMP_predictions.tsv",
                mime="text/tab-separated-values",
                use_container_width=True,
                key="dl_amp_tsv",
                disabled=amp_result_df.empty,
            )
        with col_amp_csv:
            st.download_button(
                f"📥 AMP CTD 피처.csv  ({len(amp_feat_df)}개)",
                data=_to_csv(amp_feat_df),
                file_name="ampep_AMP_features.csv",
                mime="text/csv",
                use_container_width=True,
                key="dl_amp_csv",
                disabled=amp_feat_df.empty,
            )

        # ── non-AMP ──
        lbl_namp, col_namp_tsv, col_namp_csv = st.columns([1.2, 2, 2])
        lbl_namp.markdown("🔴 **non-AMP**")
        with col_namp_tsv:
            st.download_button(
                f"📥 non-AMP 예측 결과.tsv  ({len(nonamp_result_df)}개)",
                data=_to_tsv(nonamp_result_df),
                file_name="ampep_nonAMP_predictions.tsv",
                mime="text/tab-separated-values",
                use_container_width=True,
                key="dl_namp_tsv",
                disabled=nonamp_result_df.empty,
            )
        with col_namp_csv:
            st.download_button(
                f"📥 non-AMP CTD 피처.csv  ({len(nonamp_feat_df)}개)",
                data=_to_csv(nonamp_feat_df),
                file_name="ampep_nonAMP_features.csv",
                mime="text/csv",
                use_container_width=True,
                key="dl_namp_csv",
                disabled=nonamp_feat_df.empty,
            )

