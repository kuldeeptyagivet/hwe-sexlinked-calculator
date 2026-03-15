# 🧬 HWE at Sex-Linked Loci — Interactive Calculator

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://hwe-sexlinked-calculator.streamlit.app)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)

**Companion computational tool for:**

> Tyagi KK (2026). *Approach to Hardy–Weinberg Equilibrium at Sex-Linked Loci: An Elementary Derivation with Applications to Animal Breeding.* Genetics Selection Evolution *(under review)*

---

## 📌 Overview

This browser-based interactive calculator computes allele frequency dynamics at sex-linked loci under Hardy–Weinberg equilibrium. It implements the closed-form formulae derived in the companion paper and provides practical outputs for animal breeders and genomic QC analysts.

The tool addresses a gap in existing software: **no standalone calculator existed for sex-linked loci that quantifies the approach to HWE from arbitrary initial allele frequencies in two sexes.**

---

## 🚀 Live Application

👉 **[Launch the app here](https://hwe-sexlinked-calculator.streamlit.app)**

No installation required. Runs in any modern web browser.

---

## ✨ Features

| Tab | What it does |
|-----|-------------|
| 📈 Convergence Plot | Interactive plots of female/male allele frequency trajectories with equilibrium band |
| 📊 Frequency Table | Generation-by-generation qf(t), qm(t), deviations — downloadable as CSV |
| ⏱️ Generations to Equilibrium | Closed-form t_min formula with universal reference table |
| 🔬 HWE Test Analysis | Expected HWE test failure duration for genomic QC datasets |
| 📐 Jacobsthal Numbers | Generation-by-generation Jacobsthal coefficients with convergence demonstration |

---

## 🧮 Mathematical Background

For a sex-linked locus with initial allele frequencies qf₀ (females) and qm₀ (males), random mating produces the recurrences:

```
qf(t+1) = [qf(t) + qm(t)] / 2        # females inherit X from both parents
qm(t+1) = qf(t)                        # males inherit X only from mother
```

The equilibrium frequency and closed-form solution are:

```
q̄ = (2·qf₀ + qm₀) / 3

qf(t) = q̄ + [(qf₀ − qm₀)/3] · (−½)ᵗ
qm(t) = q̄ − [2(qf₀ − qm₀)/3] · (−½)ᵗ
```

The minimum generations to reach within tolerance ε of equilibrium:

```
t_min = ⌈ log₂(2d / 3ε) ⌉    where d = |qf₀ − qm₀|
```

The coefficients in the generalized formula follow the **Jacobsthal number sequence**
(0, 1, 1, 3, 5, 11, 21, ..., OEIS A001045), with closed form Jₙ = [2ⁿ − (−1)ⁿ] / 3.

---

## 📥 Installation (Local)

To run the app locally on your machine:

```bash
# Clone the repository
git clone https://github.com/[your-username]/hwe-sexlinked-calculator.git
cd hwe-sexlinked-calculator

# Install dependencies
pip install -r requirements.txt

# Launch the app
streamlit run hwe_sexlinked_app.py
```

The app will open automatically in your browser at `http://localhost:8501`

---

## 📋 Requirements

```
streamlit>=1.28.0
pandas>=1.5.0
numpy>=1.23.0
matplotlib>=3.6.0
```

---

## 🐄 Applications in Animal Breeding

This tool is specifically designed for animal genetics and breeding contexts:

- **Crossbreeding programmes** — predict how many generations before a sex-linked locus reaches HWE after crossing two breeds with different allele frequencies
- **Genomic selection** — identify how long sex-linked SNPs will show spurious HWE deviations post-cross, preventing false removal from SNP panels
- **Breed formation** — estimate time to HWE at sex-linked loci of economic importance (coat colour, polledness, production traits)
- **Teaching** — interactive demonstration of criss-cross inheritance and Jacobsthal number patterns for undergraduate/postgraduate genetics courses

---

## 🔢 Key Results Accessible Through the Tool

1. **Convergence is rapid** — even with maximum initial disparity (d = 0.9), within 4 generations both sexes are within 5% of equilibrium
2. **Males always lag** — male deviation from equilibrium is always exactly 2× the female deviation
3. **Sex ratio does not affect trajectory** — the deterministic allele frequency path is identical regardless of female-to-male mating ratio (1:1 to 50:1)
4. **HWE failures are predictable** — sex-linked SNPs in crossbred livestock datasets fail HWE testing for at most 1–3 generations; this is expected mathematics, not genotyping error

---

## 📖 Citation

If you use this tool in your research, please cite:

```
Tyagi KK (2026). Approach to Hardy–Weinberg Equilibrium at Sex-Linked Loci:
An Elementary Derivation with Applications to Animal Breeding.
Genetics Selection Evolution (under review).
```

And for the prior pedagogical derivation:

```
Tyagi K (2021). Hardy Weinberg Law. Lecture notes, Principles of Animal and
Population Genetics AGB Unit II. Sardar Vallabhbhai Patel University of
Agriculture & Technology, Meerut, India. Retrieved from https://vepub.com
```

---

## 👤 Author

**Dr. Kuldeep Kumar Tyagi**
Associate Professor & Officer-in-Charge
Department of Animal Genetics & Breeding
College of Veterinary & Animal Sciences
Sardar Vallabhbhai Patel University of Agriculture & Technology
Meerut – 250 110, Uttar Pradesh, India

📧 drtyagivet@gmail.com
🌐 https://sites.google.com/view/learnagb

---

## 📚 Key References

- Goldberg A, Rosenberg NA (2015). Beyond 2/3 and 1/3: the complex signatures of sex-biased admixture on the X chromosome. *Genetics* 201: 263–279.
- Hosking L et al. (2004). Detection of genotyping errors by Hardy–Weinberg equilibrium testing. *European Journal of Human Genetics* 12: 395–399.
- Jennings HS (1916). The numerical results of diverse systems of breeding. *Genetics* 1: 53–89.
- Rosenberg NA (2016). Admixture models and the breeding systems of H.S. Jennings: A GENETICS connection. *Genetics* 202: 9–13.
- Crow JF, Kimura M (1970). *An Introduction to Population Genetics Theory.* Burgess, Minneapolis.

---

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

*Department of Animal Genetics & Breeding, COVAS, SVPUAT, Meerut, India*
