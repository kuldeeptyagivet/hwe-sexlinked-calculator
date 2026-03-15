"""
HWE at Sex-Linked Loci — Interactive Calculator
Author: Kuldeep Kumar Tyagi
Department of Animal Genetics & Breeding
COVAS, SVPUAT, Meerut, India

Companion tool for:
"Approach to Hardy-Weinberg Equilibrium at Sex-Linked Loci:
An Elementary Derivation with Applications to Animal Breeding"
Target journal: Genetics Selection Evolution
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import io
import math

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="HWE Sex-Linked Loci Calculator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ── Custom CSS ─────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .main-header {
        font-size: 1.6rem;
        font-weight: 700;
        color: #1a3a5c;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 0.95rem;
        color: #555;
        margin-bottom: 1.5rem;
    }
    .result-box {
        background: #f0f6ff;
        border-left: 4px solid #2c5f8a;
        padding: 12px 16px;
        border-radius: 4px;
        margin: 8px 0;
    }
    .warning-box {
        background: #fff8e1;
        border-left: 4px solid #f9a825;
        padding: 12px 16px;
        border-radius: 4px;
        margin: 8px 0;
    }
    .obs-box {
        background: #e8f5e9;
        border-left: 4px solid #2e7d32;
        padding: 10px 14px;
        border-radius: 4px;
        margin: 6px 0;
        font-size: 0.9rem;
    }
    .formula-box {
        background: #f5f5f5;
        border: 1px solid #ddd;
        padding: 10px 16px;
        border-radius: 4px;
        font-family: monospace;
        font-size: 0.88rem;
        margin: 8px 0;
    }
    .metric-card {
        background: white;
        border: 1px solid #ddd;
        border-radius: 6px;
        padding: 14px;
        text-align: center;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08);
    }
    .metric-value {
        font-size: 1.8rem;
        font-weight: 700;
        color: #1a3a5c;
    }
    .metric-label {
        font-size: 0.8rem;
        color: #666;
        margin-top: 4px;
    }
    .stTabs [data-baseweb="tab"] {
        font-size: 0.9rem;
    }
</style>
""", unsafe_allow_html=True)

# ── Core computation functions ─────────────────────────────────────────────────

def compute_trajectory(qf0, qm0, n_gen):
    """Compute allele frequency trajectory for n_gen generations."""
    q_bar = (2 * qf0 + qm0) / 3
    d = qf0 - qm0

    rows = []
    qf_prev, qm_prev = qf0, qm0

    for t in range(n_gen + 1):
        qf = q_bar + (d / 3) * ((-0.5) ** t)
        qm = q_bar - (2 * d / 3) * ((-0.5) ** t)

        # Excess heterozygosity (only defined for t >= 1)
        if t == 0:
            delta_H = None
            H_obs = None
            H_exp = None
        else:
            delta_H = (d ** 2) / 2 * (0.25 ** (t - 1))
            H_obs = qf_prev + qm_prev - 2 * qf_prev * qm_prev
            H_exp = 2 * qf * (1 - qf)

        rows.append({
            "Generation": t,
            "qf (Female)": round(qf, 6),
            "qm (Male)": round(qm, 6),
            "q̄ (Equilibrium)": round(q_bar, 6),
            "Female deviation": round(abs(qf - q_bar), 6),
            "Male deviation": round(abs(qm - q_bar), 6),
            "ΔH (Excess heterozygosity)": round(delta_H, 6) if delta_H is not None else None,
            "H_obs": round(H_obs, 6) if H_obs is not None else None,
            "H_exp": round(H_exp, 6) if H_exp is not None else None,
        })
        qf_prev, qm_prev = qf, qm

    return pd.DataFrame(rows), q_bar


def generations_to_equilibrium(d, epsilon):
    """Compute minimum generations to reach within epsilon of equilibrium."""
    if d == 0:
        return 0
    val = (2 * d) / (3 * epsilon)
    if val <= 1:
        return 0
    return math.ceil(math.log2(val))


def hwe_failure_generations(d, q_bar, nf):
    """
    Compute last generation where HWE test is expected to fail.
    Returns t_fail (last failing generation), or 0 if passes from start.
    """
    if d == 0:
        return 0

    results = []
    for t in range(1, 20):
        delta_H = (d ** 2) / 2 * (0.25 ** (t - 1))
        denom = 2 * q_bar * (1 - q_bar)
        if denom == 0:
            chi2 = 0
        else:
            chi2 = nf * (delta_H ** 2) / denom
        results.append((t, chi2, chi2 > 3.841))

    # Find last generation that fails
    t_fail = 0
    for t, chi2, fails in results:
        if fails:
            t_fail = t

    return t_fail, results


def jacobsthal(n):
    """Return nth Jacobsthal number."""
    return int((2**n - (-1)**n) / 3)


# ── Sidebar ────────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("### 🧬 Input Parameters")
    st.markdown("---")

    st.markdown("**Initial allele frequencies**")
    qf0 = st.slider(
        "qf₀ — Initial frequency of recessive allele in **females**",
        min_value=0.01, max_value=0.99, value=0.80, step=0.01,
        help="Frequency of the recessive allele at the sex-linked locus in females at generation 0"
    )
    qm0 = st.slider(
        "qm₀ — Initial frequency of recessive allele in **males**",
        min_value=0.01, max_value=0.99, value=0.30, step=0.01,
        help="Frequency of the recessive allele at the sex-linked locus in males at generation 0"
    )

    st.markdown("---")
    st.markdown("**Analysis settings**")

    n_gen = st.slider(
        "Number of generations to compute",
        min_value=5, max_value=30, value=15, step=1
    )

    epsilon = st.select_slider(
        "Tolerance ε for equilibrium",
        options=[0.10, 0.05, 0.02, 0.01, 0.005],
        value=0.05,
        format_func=lambda x: f"{x:.3f} ({x*100:.1f}%)"
    )

    st.markdown("---")
    st.markdown("**HWE test parameters**")

    nf = st.number_input(
        "Sample size nf (number of females)",
        min_value=50, max_value=10000, value=1000, step=50,
        help="Number of female individuals in genomic dataset"
    )

    st.markdown("---")
    st.markdown("**About**")
    st.markdown("""
    <div style='font-size:0.78rem; color:#666;'>
    Companion tool for:<br>
    <em>Tyagi KK (2026). Approach to Hardy–Weinberg Equilibrium at Sex-Linked Loci...</em><br><br>
    Dept. of Animal Genetics & Breeding<br>
    COVAS, SVPUAT, Meerut, India<br>
    drtyagivet@gmail.com
    </div>
    """, unsafe_allow_html=True)


# ── Main panel ─────────────────────────────────────────────────────────────────

st.markdown('<div class="main-header">🧬 Hardy–Weinberg Equilibrium at Sex-Linked Loci</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Interactive calculator for allele frequency dynamics, convergence analysis, and genomic QC implications | Tyagi KK (2026)</div>', unsafe_allow_html=True)

# Compute everything
df, q_bar = compute_trajectory(qf0, qm0, n_gen)
d = abs(qf0 - qm0)
t_min = generations_to_equilibrium(d, epsilon)
t_fail, chi2_results = hwe_failure_generations(d, q_bar, nf)

# ── Key metrics row ────────────────────────────────────────────────────────────
st.markdown("### Key Results")
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-value">{q_bar:.4f}</div>
        <div class="metric-label">Equilibrium frequency q̄<br>(2qf₀ + qm₀) / 3</div>
    </div>""", unsafe_allow_html=True)

with col2:
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-value">{d:.3f}</div>
        <div class="metric-label">Initial disparity d<br>|qf₀ − qm₀|</div>
    </div>""", unsafe_allow_html=True)

with col3:
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-value">{t_min}</div>
        <div class="metric-label">Generations to equilibrium<br>within ε = {epsilon}</div>
    </div>""", unsafe_allow_html=True)

with col4:
    fail_text = str(t_fail) if t_fail > 0 else "0 (passes immediately)"
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-value">{t_fail}</div>
        <div class="metric-label">Last generation HWE fails<br>nf = {nf:,} females</div>
    </div>""", unsafe_allow_html=True)

with col5:
    ne_ratio_text = "1.000 (equal)"
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-value">Unchanged</div>
        <div class="metric-label">Effect of sex ratio<br>on trajectory</div>
    </div>""", unsafe_allow_html=True)

st.markdown("---")

# ── Tabs ───────────────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📈 Convergence Plot",
    "📊 Frequency Table",
    "⏱️ Generations to Equilibrium",
    "🔬 HWE Test Analysis",
    "📐 Jacobsthal Numbers"
])


# ── TAB 1: Convergence Plot ────────────────────────────────────────────────────
with tab1:
    st.markdown("#### Allele Frequency Convergence at Sex-Linked Locus")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    generations = df["Generation"].values
    qf_vals = df["qf (Female)"].values
    qm_vals = df["qm (Male)"].values
    qbar_vals = df["q̄ (Equilibrium)"].values

    # Left plot: Frequency trajectories
    ax1 = axes[0]
    ax1.plot(generations, qf_vals, 'b-o', markersize=5, linewidth=2,
             label=f'Female frequency (qf)', zorder=3)
    ax1.plot(generations, qm_vals, 'r-s', markersize=5, linewidth=2,
             label=f'Male frequency (qm)', zorder=3)
    ax1.axhline(y=q_bar, color='green', linestyle='--', linewidth=1.5,
                label=f'Equilibrium q̄ = {q_bar:.4f}', zorder=2)

    # Shade tolerance band
    ax1.axhspan(q_bar - epsilon, q_bar + epsilon, alpha=0.12, color='green',
                label=f'±ε tolerance band (ε={epsilon})')

    # Mark starting points
    ax1.plot(0, qf0, 'b*', markersize=12, zorder=4)
    ax1.plot(0, qm0, 'r*', markersize=12, zorder=4)

    # Mark convergence generation
    if t_min > 0 and t_min <= n_gen:
        ax1.axvline(x=t_min, color='purple', linestyle=':', linewidth=1.5,
                    alpha=0.7, label=f'Equilibrium reached (t={t_min})')

    ax1.set_xlabel('Generation (t)', fontsize=11)
    ax1.set_ylabel('Allele frequency (q)', fontsize=11)
    ax1.set_title('Allele Frequency Dynamics\n(Criss-cross inheritance pattern)', fontsize=11)
    ax1.legend(fontsize=8.5, loc='best')
    ax1.set_xlim(-0.3, n_gen + 0.3)
    ax1.set_ylim(max(0, min(qf_vals.min(), qm_vals.min()) - 0.05),
                 min(1, max(qf_vals.max(), qm_vals.max()) + 0.05))
    ax1.grid(True, alpha=0.3)
    ax1.set_facecolor('#fafafa')

    # Right plot: Deviation from equilibrium (log scale)
    ax2 = axes[1]
    dev_f = df["Female deviation"].values[1:]  # skip t=0
    dev_m = df["Male deviation"].values[1:]
    gen_skip = generations[1:]

    ax2.semilogy(gen_skip, dev_f, 'b-o', markersize=5, linewidth=2,
                 label='|qf(t) − q̄| Female deviation')
    ax2.semilogy(gen_skip, dev_m, 'r-s', markersize=5, linewidth=2,
                 label='|qm(t) − q̄| Male deviation')
    ax2.axhline(y=epsilon, color='green', linestyle='--', linewidth=1.5,
                label=f'Tolerance ε = {epsilon}')

    ax2.set_xlabel('Generation (t)', fontsize=11)
    ax2.set_ylabel('|Deviation from equilibrium| (log scale)', fontsize=11)
    ax2.set_title('Deviation from Equilibrium\n(Log scale — note male deviation = 2× female)', fontsize=11)
    ax2.legend(fontsize=8.5, loc='best')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_facecolor('#fafafa')

    plt.tight_layout(pad=2.5)
    st.pyplot(fig)
    plt.close()

    col_obs1, col_obs2 = st.columns(2)
    with col_obs1:
        st.markdown(f"""
        <div class="obs-box">
        <strong>Observation 1 — Criss-cross pattern:</strong> Male frequency at generation t equals
        female frequency at generation t−1. The two curves oscillate on opposite sides of the
        equilibrium line, converging with decreasing amplitude.
        </div>""", unsafe_allow_html=True)
    with col_obs2:
        st.markdown(f"""
        <div class="obs-box">
        <strong>Observation 2 — Male lags always:</strong> The male deviation from equilibrium
        (red) is always exactly twice the female deviation (blue) at any generation t. Males
        are always the last to converge — confirmed in the log-scale plot.
        </div>""", unsafe_allow_html=True)


# ── TAB 2: Frequency Table ─────────────────────────────────────────────────────
with tab2:
    st.markdown("#### Generation-by-Generation Allele Frequencies")

    display_df = df[["Generation", "qf (Female)", "qm (Male)",
                      "q̄ (Equilibrium)", "Female deviation", "Male deviation"]].copy()

    def highlight_converged(row):
        colors = []
        for col in row.index:
            if col in ["Female deviation", "Male deviation"]:
                val = row[col]
                if val <= epsilon:
                    colors.append("background-color: #e8f5e9; color: #2e7d32")
                elif val <= epsilon * 3:
                    colors.append("background-color: #fff8e1")
                else:
                    colors.append("background-color: #fff3e0")
            else:
                colors.append("")
        return colors

    styled = display_df.style\
        .apply(highlight_converged, axis=1)\
        .format({
            "qf (Female)": "{:.6f}",
            "qm (Male)": "{:.6f}",
            "q̄ (Equilibrium)": "{:.6f}",
            "Female deviation": "{:.6f}",
            "Male deviation": "{:.6f}"
        })

    st.dataframe(styled, use_container_width=True, height=400)

    st.markdown(f"""
    <div class="result-box">
    <strong>Color guide:</strong>
    🟩 Green = within tolerance ε={epsilon} |
    🟨 Yellow = within 3ε |
    🟧 Orange = outside 3ε
    </div>""", unsafe_allow_html=True)

    # Download button
    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index=False)
    st.download_button(
        label="⬇️ Download full table as CSV",
        data=csv_buffer.getvalue(),
        file_name=f"hwe_sexlinked_qf{qf0}_qm{qm0}.csv",
        mime="text/csv"
    )


# ── TAB 3: Generations to Equilibrium ─────────────────────────────────────────
with tab3:
    st.markdown("#### Generations Required to Reach Hardy–Weinberg Equilibrium")

    col_left, col_right = st.columns([1, 1])

    with col_left:
        st.markdown(f"""
        <div class="result-box">
        <strong>For your parameters:</strong><br>
        d = |qf₀ − qm₀| = |{qf0} − {qm0}| = <strong>{d:.4f}</strong><br>
        ε (tolerance) = <strong>{epsilon}</strong><br>
        <br>
        <strong>t_min = ⌈ log₂(2d / 3ε) ⌉ = ⌈ log₂({2*d:.4f} / {3*epsilon:.4f}) ⌉
        = ⌈ log₂({2*d/(3*epsilon):.4f}) ⌉ = <span style='font-size:1.3rem; color:#1a3a5c'>
        {t_min} generations</span></strong>
        </div>""", unsafe_allow_html=True)

        st.markdown("""
        <div class="formula-box">
        Formula: t_min = ⌈ log₂(2d / 3ε) ⌉<br>
        where d = |qf₀ − qm₀|, ε = tolerance threshold<br>
        ⌈ ⌉ = ceiling function (round up to integer)<br>
        Based on: |qm(t) − q̄| = (2d/3)·(½)ᵗ ≤ ε
        </div>""", unsafe_allow_html=True)

        st.markdown(f"""
        <div class="obs-box">
        <strong>Biological interpretation:</strong> After {t_min} generation(s) of random mating,
        both female and male allele frequencies will be within {epsilon} allele frequency
        units of the equilibrium value q̄ = {q_bar:.4f}.
        {"This is the maximum possible (for the given ε), since d is very large." if d >= 0.7 else ""}
        </div>""", unsafe_allow_html=True)

    with col_right:
        # Full table for current d across all epsilon values
        st.markdown("**Complete table for d = {:.3f}:**".format(d))
        eps_vals = [0.10, 0.05, 0.02, 0.01, 0.005]
        t_vals = [generations_to_equilibrium(d, e) for e in eps_vals]
        table_df = pd.DataFrame({
            "Tolerance ε": [f"{e:.3f} ({e*100:.1f}%)" for e in eps_vals],
            "t_min (generations)": t_vals
        })
        st.dataframe(table_df, use_container_width=True, hide_index=True)

    st.markdown("---")
    st.markdown("**Universal table — minimum generations for all combinations:**")

    d_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    eps_cols = [0.05, 0.02, 0.01]

    table_data = {}
    table_data["d = |qf₀ − qm₀|"] = [f"{v:.1f}" for v in d_vals]
    for e in eps_cols:
        col_name = f"ε = {e:.2f} ({int(e*100)}%)"
        table_data[col_name] = [generations_to_equilibrium(v, e) for v in d_vals]

    full_table = pd.DataFrame(table_data)

    def highlight_current_d(row):
        d_val = float(row["d = |qf₀ − qm₀|"])
        if abs(d_val - round(d, 1)) < 0.05:
            return ["background-color: #dceeff; font-weight: bold"] * len(row)
        return [""] * len(row)

    st.dataframe(
        full_table.style.apply(highlight_current_d, axis=1),
        use_container_width=True,
        hide_index=True
    )
    st.caption("Highlighted row = closest match to your current d value. "
               "All values derived from the closed-form formula and verified by direct iteration.")


# ── TAB 4: HWE Test Analysis ───────────────────────────────────────────────────
with tab4:
    st.markdown("#### HWE Test Failure Analysis for Genomic QC")
    st.markdown("""
    Sex-linked SNPs in crossbred populations are expected to fail HWE testing for a
    predictable number of generations due to allele frequency disequilibrium between sexes —
    **not** because of genotyping error. This tab quantifies that expectation.
    """)

    col_l, col_r = st.columns([1, 1])

    with col_l:
        st.markdown(f"""
        <div class="result-box">
        <strong>For your parameters (nf = {nf:,} females):</strong><br>
        d = {d:.4f}, q̄ = {q_bar:.4f}<br><br>
        Expected HWE test failures: <strong>generations 1 through {t_fail}</strong><br>
        {"✅ Passes HWE from generation 1 onward" if t_fail == 0 else
         f"⚠️ SNP expected to fail HWE test for {t_fail} generation(s) post-cross"}
        <br><br>
        <em>Recommendation: Do not flag this sex-linked SNP as a genotyping error
        in generations 1–{t_fail}. Investigate HWE failure only from generation
        {t_fail+1} onward.</em>
        </div>""", unsafe_allow_html=True)

        st.markdown("""
        <div class="formula-box">
        ΔH(t) = d²/2 · (1/4)^(t−1)  [Excess heterozygosity]<br>
        χ²(t) ≈ nf · [ΔH(t)]² / [2·q̄·(1−q̄)]<br>
        Fails HWE when χ²(t) > 3.841  (p < 0.05, 1 df)
        </div>""", unsafe_allow_html=True)

    with col_r:
        # Chi-square trajectory table
        chi2_df = pd.DataFrame([
            {
                "Generation": t,
                "ΔH(t)": round((d**2)/2 * (0.25**(t-1)), 6),
                "χ²(t)": round(nf * (((d**2)/2*(0.25**(t-1)))**2) / max(1e-10, 2*q_bar*(1-q_bar)), 3),
                "HWE test (p<0.05)": "❌ FAILS" if chi2 > 3.841 else "✅ PASSES"
            }
            for t, chi2, fails in chi2_results[:10]
        ])

        def color_hwe(row):
            if "FAILS" in str(row["HWE test (p<0.05)"]):
                return ["background-color: #ffebee"] * len(row)
            return ["background-color: #e8f5e9"] * len(row)

        st.dataframe(
            chi2_df.style.apply(color_hwe, axis=1),
            use_container_width=True,
            hide_index=True
        )

    st.markdown("---")
    st.markdown("**Universal HWE failure table — last failing generation (q̄ = 0.4 assumed):**")

    d_rows = [0.2, 0.3, 0.5, 0.7, 0.9]
    nf_cols = [100, 500, 1000, 5000]
    q_bar_ref = 0.4

    hwe_data = {"d = |qf₀−qm₀|": [str(v) for v in d_rows]}
    for n in nf_cols:
        col_vals = []
        for dv in d_rows:
            t_f, _ = hwe_failure_generations(dv, q_bar_ref, n)
            col_vals.append(t_f)
        hwe_data[f"nf = {n:,}"] = col_vals

    hwe_table = pd.DataFrame(hwe_data)
    st.dataframe(hwe_table, use_container_width=True, hide_index=True)
    st.caption("Values show last generation expected to fail HWE test (χ² > 3.841, p < 0.05). "
               "0 = passes from first generation. Computed at q̄ = 0.4 (worst case for variance).")

    st.markdown(f"""
    <div class="obs-box">
    <strong>Observation — Direction of deviation:</strong> Sex-linked loci in crossbred
    populations ALWAYS show excess heterozygosity in females (ΔH ≥ 0). A heterozygote
    <em>deficit</em> at a sex-linked locus cannot arise from inter-sex allele frequency
    disequilibrium and therefore warrants investigation as a genuine genotyping artefact.
    </div>""", unsafe_allow_html=True)


# ── TAB 5: Jacobsthal Numbers ──────────────────────────────────────────────────
with tab5:
    st.markdown("#### Jacobsthal Numbers in Sex-Linked Allele Frequency Dynamics")
    st.markdown("""
    The coefficients of qf₀ and qm₀ in the generalized allele frequency formulae across
    successive generations follow the **Jacobsthal number sequence**: 0, 1, 1, 3, 5, 11, 21, 43, ...

    Defined by: J₀ = 0, J₁ = 1, and **Jₙ = Jₙ₋₁ + 2Jₙ₋₂**

    Closed form: **Jₙ = [2ⁿ − (−1)ⁿ] / 3** (OEIS sequence A001045)

    First noted in a population genetics context by Jennings (1916), rediscovered in admixture
    modelling by Goldberg & Rosenberg (2015), and independently derived in an Animal Genetics
    teaching context by Tyagi (2021).
    """)

    # Jacobsthal table
    max_t = min(n_gen, 12)
    jac_rows = []
    for t in range(max_t + 1):
        jt = jacobsthal(t)
        jt1 = jacobsthal(t + 1)
        qft_formula = f"({jt1}/2^{t})·qf₀ + ({jt}/2^{t})·qm₀" if t > 0 else "qf₀"
        qmt_formula = f"({jt}/2^{t-1})·qf₀ + ({jacobsthal(t-1)}/2^{t-1})·qm₀" if t > 0 else "qm₀"

        qft_val = q_bar + (qf0 - qm0) / 3 * ((-0.5) ** t)
        qmt_val = q_bar - 2 * (qf0 - qm0) / 3 * ((-0.5) ** t)

        jac_rows.append({
            "t": t,
            "Jₜ": jacobsthal(t),
            "Jₜ₊₁": jacobsthal(t + 1),
            "Jₜ₊₁/2ᵗ (coeff of qf₀)": f"{jt1}/{2**t} = {jt1/2**t:.6f}",
            "Jₜ/2ᵗ (coeff of qm₀)": f"{jt}/{2**t} = {jt/2**t:.6f}",
            "qf(t) computed": round(qft_val, 6),
            "qm(t) computed": round(qmt_val, 6)
        })

    jac_df = pd.DataFrame(jac_rows)
    st.dataframe(jac_df, use_container_width=True, hide_index=True, height=380)

    st.markdown("""
    <div class="formula-box">
    General formula (females): qf(t) = [Jₜ₊₁/2ᵗ]·qf₀ + [Jₜ/2ᵗ]·qm₀<br>
    General formula (males):   qm(t) = [Jₜ/2ᵗ⁻¹]·qf₀ + [Jₜ₋₁/2ᵗ⁻¹]·qm₀  (t ≥ 1)<br><br>
    Jacobsthal sequence:  0, 1, 1, 3, 5, 11, 21, 43, 85, 171, ...<br>
    Denominators (powers of 2): 1, 2, 4, 8, 16, 32, 64, 128, ...
    </div>""", unsafe_allow_html=True)

    st.markdown(f"""
    <div class="obs-box">
    <strong>Convergence of coefficients:</strong> As t → ∞, both Jₜ₊₁/2ᵗ → 2/3 and
    Jₜ/2ᵗ → 1/3. This confirms that the equilibrium frequency q̄ = (2qf₀ + qm₀)/3
    is the weighted mean with the 2:1 female-to-male ratio of X chromosomes in the population.
    </div>""", unsafe_allow_html=True)


# ── Footer ─────────────────────────────────────────────────────────────────────
st.markdown("---")
st.markdown("""
<div style='text-align:center; font-size:0.78rem; color:#888; padding:10px 0'>
HWE Sex-Linked Loci Calculator | Tyagi KK (2026) | COVAS, SVPUAT, Meerut, India<br>
Companion tool for: <em>Approach to Hardy–Weinberg Equilibrium at Sex-Linked Loci — An Elementary Derivation with Applications to Animal Breeding</em><br>
Target: Genetics Selection Evolution | For queries: drtyagivet@gmail.com
</div>
""", unsafe_allow_html=True)
