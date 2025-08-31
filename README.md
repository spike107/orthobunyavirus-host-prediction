# Orthobunyavirus Host-Prediction Pipeline

> MSc Bioinformatics project — reproducible scripts for building virus host labels, extracting genomic/protein features, and training machine-learning models (Random Forest) under cross-validation and leave-one-virus-out evaluation.

---
⚠️ Disclaimer  
This repository is a **research prototype** accompanying an MSc dissertation.  
It is not production-ready software and should not be used for clinical or public health decision-making.  
Results may vary depending on dataset versions and parameter choices.

## 1. Overview

This repository contains a complete implementation of the pipeline described in the dissertation:
- **Phase 1 — Label table construction**: ICTV/DB-based standardisation, evidence-aware relabelling (human vs nonhuman).
- **Phase 2 — FASTA preprocessing**: parse headers, strip version suffixes, build clean IDs, map across nucleotide/CDS/protein sets.
- **Phase 3 — Feature engineering and ML**: nucleotide k-mer (k=2–6), codon usage (RSCU, CAI), protein descriptors (via iFeatureOmegaCLI).  
  Models trained under stratified 80/20 and leave-one-virus-out (LOVO) with group-aware CV parameter tuning.

Key script: `src/project_pipeline.py`.

---

## 2. Repository Layout

```
.
├─ README.md
├─ environment.yml
├─ requirements.txt
├─ LICENSE
├─ CITATION.cff
├─ .gitignore
├─ data/            # optional: place inputs here if using README-style names
├─ results/         # created automatically; all outputs collected here
└─ src/
   └─ project_pipeline.py
```

> Note: the `data/` directory is optional. Users may either place input files in the repository root (legacy names) or under `data/` (README-style names).  
> The `results/` directory is created automatically during runtime.

---


## 3. Installation

### Conda (recommended)
```bash
conda env create -f environment.yml
conda activate orthobunya
```

### pip
```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### External tools
- [iFeatureOmegaCLI](https://github.com/Superzchen/iFeatureOmega) for protein descriptors.  
  This tool is **optional**: if not installed or if `--protein` is unset, protein features are skipped gracefully.  
  The pipeline falls back to available feature sets (priority: **RSCU → k-mer → protein**).

---

## 4. Data acquisition

Viral sequences were obtained from the **NCBI Virus portal** (Nucleotide, CDS, and Protein databases), retrieving **all records classified under the *Orthobunyavirus* genus** (NCBI Taxon ID: 11572).  
No filtering was applied at the download stage; instead, completeness checks and length thresholds (e.g. CDS ≥600 nt for RSCU/CAI) were enforced downstream in the pipeline.

Raw sequence data are **not included** in this repository due to size and licensing constraints.  
Users should download the required FASTA and reference tables (see *Inputs*) and place them either in the repository root (legacy names) or under the `data/` directory (README-style names).

---

## 5. Inputs

The pipeline expects several input FASTA and reference files.  
Two naming conventions are supported:

### A. Legacy names (place directly in the repository root)
- `sequences_human.fasta`  
- `sequences_nothuman.fasta`  
- `sequences_CDS.fasta`  
- `sequences_protein.fasta`  
- `human_HK_CDS.cleaned.fasta`  
- `Orthobunyavirus_taxonomy.xlsx`  
- `ICTV_2024_Taxa Renamed or Abolished.csv`  
- `ICTV_2024_MSL.csv`  
- `ICTV Orthobunyavirus.csv`  
- `human_virus_DB_species_clean.txt`  

### B. README-style names (place in `data/` subdirectory)
- `data/nucleotide.fasta`  
- `data/nonhuman.fasta`  
- `data/cds.fasta`  
- `data/protein.fasta`  
- `data/human_HK_CDS.cleaned.fasta`  
- Other reference tables may also be placed under `data/`

---

## 6. Usage

### Windows example
```powershell
# If 'conda activate' is not available, call python.exe inside your conda env:
C:/Users/<user>/.conda/envs/<env_name>/python.exe src/project_pipeline.py --help

# Minimal run (auto-detect inputs; outputs to results_cli/)
C:/Users/<user>/.conda/envs/<env_name>/python.exe src/project_pipeline.py --outdir results_cli
```

### Linux / macOS example
```bash
conda activate orthobunya
python src/project_pipeline.py --outdir results_cli
```

Notes:
- If inputs are not provided via CLI, the script automatically detects files in the **repo root (legacy names)**, falling back to **`data/` (README names)**.  
- All new outputs are collected under `results/` (or `--outdir`). Logs are saved there as `pipeline_log.txt`.

---

## 7. Results

- All outputs generated at the repository root are automatically moved to `results/` (or the directory passed with `--outdir`).  
- Subfolders include:  
  - `tables/`: metadata, train-test and LOVO results, feature-importance summaries.  
  - `figs/`: performance plots, LOVO boxplots, UMAPs, clustermaps, feature-importance plots.  
  - `logs/`: runtime logs.  

---

## 8. Computational environment and reproducibility

Analyses were run in:
- Python ≥3.8  
- Biopython ≥1.8  
- pandas ≥1.5  
- scikit-learn ≥1.2  
- matplotlib ≥3.7  
- seaborn ≥0.12  
- umap-learn ≥0.5  
- iFeatureOmegaCLI (optional external tool)  

Randomised procedures used fixed seeds (42):
- Random Forest (`random_state=42`)  
- Stratified 80/20 split (`random_state=42`)  
- Group-aware 5-fold CV with shuffling (`StratifiedGroupKFold(..., random_state=42)`)  

Intermediate outputs, logs, and provenance files are retained for reproducibility.

---

## 9. Intersection fallback

- When all three feature sets are available: **3-way intersection**.  
- When only two are available: **2-way intersection**.  
- If intersection is empty: fall back to the richest available single set (priority: **RSCU → k-mer → protein**).  
- This ensures the ML evaluation always runs, even without protein features.

---

## 10. Citation

```bibtex
@software{jin_orthobunyavirus_pipeline_2025,
  author  = {Jin, Fengshi},
  title   = {Orthobunyavirus Host-Prediction Pipeline},
  year    = {2025},
  url     = {https://github.com/spike107/orthobunyavirus-host-prediction/releases/tag/v1.0.1},
  version = {v1.0.1}
}
```

---

## 11. License

MIT License (see `LICENSE`).  
Input data are **not distributed** — users must obtain FASTA and reference tables from ICTV, Virus-Host DB, NCBI, etc.

---

## 12. Contact

Please open GitHub issues or pull requests.
