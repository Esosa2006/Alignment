import streamlit as st

st.set_page_config(
    page_title="Needleman-Wunsch Aligner",
    page_icon="🔗",
    layout="centered"
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;500&display=swap');

html, body, [class*="css"] { font-family: 'DM Sans', sans-serif; }
.stApp { background: #0f0a18; color: #e8e2f2; }
h1,h2,h3 { font-family: 'Space Mono', monospace !important; }

.main-title {
    font-family: 'Space Mono', monospace;
    font-size: 2rem; font-weight: 700;
    color: #c77df4; letter-spacing: -1px; margin-bottom: 0.2rem;
}
.sub-title { font-size: 0.95rem; color: #555; margin-bottom: 2rem; font-weight: 300; }

.result-section {
    background: #0c0914; border: 1px solid #1e1530;
    border-radius: 4px; padding: 1.25rem 1.5rem; margin-bottom: 1rem;
}
.section-label {
    font-family: 'Space Mono', monospace; font-size: 0.7rem;
    color: #444; letter-spacing: 2px; text-transform: uppercase; margin-bottom: 0.6rem;
}
.alignment-block {
    font-family: 'Space Mono', monospace; font-size: 1.1rem;
    line-height: 2; word-break: break-all;
}
.seq-purple { color: #c77df4; }
.seq-teal   { color: #7df4c8; }
.bar-char   { color: #555; }

.stat-grid {
    display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin-bottom: 1rem;
}
.stat-card {
    background: #0c0914; border: 1px solid #1e1530;
    border-radius: 4px; padding: 0.9rem 1rem; text-align: center;
}
.stat-num {
    font-family: 'Space Mono', monospace; font-size: 1.6rem;
    font-weight: 700; color: #c77df4;
}
.stat-label {
    font-family: 'Space Mono', monospace; font-size: 0.65rem;
    color: #444; letter-spacing: 1.5px; text-transform: uppercase; margin-top: 3px;
}
.score-positive { color: #7df4c8 !important; }
.score-negative { color: #f47d7d !important; }

.params-row {
    display: flex; gap: 8px; margin-bottom: 1rem; flex-wrap: wrap;
}
.param-pill {
    background: #1a1025; border: 1px solid #2a1a40;
    border-radius: 3px; padding: 4px 10px;
    font-family: 'Space Mono', monospace; font-size: 0.75rem; color: #7a5a9a;
}

.stButton > button {
    background: #c77df4 !important; color: #0f0a18 !important;
    border: none !important; border-radius: 4px !important;
    font-family: 'Space Mono', monospace !important; font-weight: 700 !important;
    font-size: 0.85rem !important; letter-spacing: 1px !important;
    padding: 0.6rem 2rem !important; width: 100% !important; margin-top: 0.5rem !important;
}
.stButton > button:hover { background: #d49af6 !important; }

.stTextInput input {
    background: #0c0914 !important; border-color: #1e1530 !important;
    color: #e8e2f2 !important; border-radius: 4px !important;
    font-family: 'Space Mono', monospace !important; font-size: 0.95rem !important;
}
.stTextInput label {
    font-family: 'Space Mono', monospace !important; font-size: 0.75rem !important;
    letter-spacing: 1.5px !important; text-transform: uppercase !important; color: #555 !important;
}
.stNumberInput input {
    background: #0c0914 !important; border-color: #1e1530 !important;
    color: #e8e2f2 !important; border-radius: 4px !important;
    font-family: 'Space Mono', monospace !important;
}
.stNumberInput label {
    font-family: 'Space Mono', monospace !important; font-size: 0.7rem !important;
    letter-spacing: 1.5px !important; text-transform: uppercase !important; color: #555 !important;
}
.error-box {
    background: #1a0808; border: 1px solid #5a1a1a;
    border-radius: 4px; padding: 1rem 1.5rem;
    font-family: 'Space Mono', monospace; font-size: 0.85rem; color: #f57a7a;
}
.legend-row {
    display: flex; gap: 16px; margin-top: 0.5rem;
}
.legend-item {
    font-family: 'Space Mono', monospace; font-size: 0.8rem; color: #555;
}
</style>
""", unsafe_allow_html=True)

def run_alignment(seq_a, seq_b, match_score, mismatch_score, gap_score):
    rows, cols = len(seq_a), len(seq_b)
    matrix = [[i*gap_score if j==0 else j*gap_score if i==0 else 0 for j in range(cols+1)] for i in range(rows+1)]
    for i in range(1, rows+1):
        for j in range(1, cols+1):
            diag = matrix[i-1][j-1] + (match_score if seq_a[i-1]==seq_b[j-1] else mismatch_score)
            matrix[i][j] = max(diag, matrix[i-1][j]+gap_score, matrix[i][j-1]+gap_score)

    aligned_a, aligned_b = [], []
    i, j = rows, cols
    while i > 0 or j > 0:
        if i>0 and j>0:
            diag = matrix[i-1][j-1] + (match_score if seq_a[i-1]==seq_b[j-1] else mismatch_score)
            if matrix[i][j] == diag:
                aligned_a.append(seq_a[i-1]); aligned_b.append(seq_b[j-1]); i-=1; j-=1; continue
        if i>0 and matrix[i][j] == matrix[i-1][j]+gap_score:
            aligned_a.append(seq_a[i-1]); aligned_b.append("-"); i-=1
        else:
            aligned_a.append("-"); aligned_b.append(seq_b[j-1]); j-=1

    result_a = "".join(reversed(aligned_a))
    result_b = "".join(reversed(aligned_b))
    match_bar = "".join("|" if x==y and x!="-" else " " if "-" in (x,y) else ":" for x,y in zip(result_a,result_b))
    return result_a, result_b, match_bar, matrix[rows][cols]

st.markdown('<div class="main-title">🔗 Needleman-Wunsch</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Global sequence alignment using dynamic programming</div>', unsafe_allow_html=True)

col1, col2 = st.columns(2)
with col1:
    seq1 = st.text_input("Sequence 1", placeholder="e.g. GATTACA")
with col2:
    seq2 = st.text_input("Sequence 2", placeholder="e.g. GCATGCU")

st.markdown("**Scoring parameters**")
pc1, pc2, pc3 = st.columns(3)
with pc1:
    match_score = st.number_input("Match reward", value=1, step=1)
with pc2:
    mismatch_score = st.number_input("Mismatch penalty", value=-1, step=1)
with pc3:
    gap_score = st.number_input("Gap penalty", value=-2, step=1)

if st.button("ALIGN"):
    s1 = seq1.strip().upper().replace(" ","")
    s2 = seq2.strip().upper().replace(" ","")

    if not s1 or not s2:
        st.markdown('<div class="error-box">⚠ Please enter both sequences.</div>', unsafe_allow_html=True)
    else:
        result_a, result_b, match_bar, score = run_alignment(s1, s2, match_score, mismatch_score, gap_score)

        matches    = match_bar.count("|")
        mismatches = match_bar.count(":")
        gaps       = result_a.count("-") + result_b.count("-")
        length     = len(result_a)
        score_cls  = "score-positive" if score >= 0 else "score-negative"

        st.markdown(f"""
        <div class="stat-grid">
            <div class="stat-card">
                <div class="stat-num {score_cls}">{score:+}</div>
                <div class="stat-label">Score</div>
            </div>
            <div class="stat-card">
                <div class="stat-num" style="color:#7df4c8">{matches}</div>
                <div class="stat-label">Matches</div>
            </div>
            <div class="stat-card">
                <div class="stat-num" style="color:#f4c47d">{mismatches}</div>
                <div class="stat-label">Mismatches</div>
            </div>
            <div class="stat-card">
                <div class="stat-num" style="color:#f47d7d">{gaps}</div>
                <div class="stat-label">Gaps</div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        def colour_seq(seq, col):
            return "".join(
                f'<span style="color:#f47d7d">-</span>' if ch == "-"
                else f'<span style="{col}">{ch}</span>'
                for ch in seq
            )

        def colour_bar(bar):
            out = ""
            for ch in bar:
                if ch == "|":   out += f'<span style="color:#7df4c8">|</span>'
                elif ch == ":": out += f'<span style="color:#f4c47d">:</span>'
                else:           out += '<span style="color:#2a2a2a">·</span>'
            return out

        st.markdown(f"""
        <div class="result-section">
            <div class="section-label">Optimal alignment (length {length})</div>
            <div class="alignment-block">{colour_seq(result_a, "color:#c77df4")}</div>
            <div class="alignment-block">{colour_bar(match_bar)}</div>
            <div class="alignment-block">{colour_seq(result_b, "color:#7df4c8")}</div>
            <div class="legend-row" style="margin-top:10px">
                <span class="legend-item"><span style="color:#7df4c8">|</span> match</span>
                <span class="legend-item"><span style="color:#f4c47d">:</span> mismatch</span>
                <span class="legend-item"><span style="color:#f47d7d">-</span> gap</span>
            </div>
        </div>
        <div class="result-section">
            <div class="section-label">Scoring used</div>
            <div class="params-row">
                <span class="param-pill">match +{match_score}</span>
                <span class="param-pill">mismatch {mismatch_score}</span>
                <span class="param-pill">gap {gap_score}</span>
                <span class="param-pill">identity {matches/length*100:.1f}%</span>
            </div>
        </div>
        """, unsafe_allow_html=True)
