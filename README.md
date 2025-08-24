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

> Note: the `config/`, `data/`, and `results/` directories are **not included** by default (they only contained placeholders).  
> They will be created locally when running the pipeline.

```
.
├─ README.md
├─ environment.yml
├─ requirements.txt
├─ LICENSE
├─ CITATION.cff
├─ .gitignore
└─ src/
   └─ project_pipeline.py
```

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
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate
pip install -r requirements.txt
```

### External tools
- iFeatureOmegaCLI (for protein descriptors).
- UMAP-learn for dimensionality reduction plots.

---

## 4. Data acquisition

Viral sequences were obtained from the **NCBI Virus portal** (Nucleotide, CDS, and Protein databases), retrieving **all records classified under the *Orthobunyavirus* genus** (NCBI Taxon ID: 11572). No filtering was applied at the download stage; instead, completeness checks and length thresholds (e.g. CDS ≥600 nt for RSCU/CAI) were enforced downstream in the pipeline.

Raw sequence data are **not included** in this repository due to size and licensing constraints. Users should download the required FASTA files from NCBI and place them under a local `data/` directory (see Inputs).

---

## 5. Inputs

Expected under a local `data/` directory (you may choose any location; paths can be overridden via command-line arguments):
- `nucleotide.fasta` — nucleotide sequences with pipe-delimited headers.
- `nonhuman.fasta` — non-human sequences if kept separate (optional; otherwise inferred).
- `cds.fasta` — CDS records; RSCU/CAI retain `complete_cds`; keep `partial_cds` ≥ 600 nt; filter `unknown` by length.
- `protein.fasta` — protein sequences (prefer those from complete genomes).
- `ictv_species_list.csv` — ICTV species list.
- `ictv_renaming_table.csv` — ICTV old↔new names.
- `virus_hostdb_human.csv` — Virus-Host DB list of human viruses.
- `manual_mapping.csv` — curated mapping.

All input file paths are provided via command-line arguments. Data files **do not** need to sit in the same folder as the script, but should be referenced consistently (e.g., `--cds data/cds.fasta`). The recommended setup is to keep all inputs in a project-level `data/` directory.

---

## 6. Usage

### Default run (with files placed under `data/`)
```bash
python src/project_pipeline.py
```
This uses the following defaults:
- `--nucleotide data/nucleotide.fasta`
- `--nonhuman  data/nonhuman.fasta`
- `--cds       data/cds.fasta`
- `--protein   data/protein.fasta`
- `--reference data/human_HK_CDS.cleaned.fasta`
- `--outdir    results/`

### Custom run (specify file paths explicitly)
```bash
python src/project_pipeline.py \
  --nucleotide mydata/human.fasta \
  --nonhuman  mydata/nonhuman.fasta \
  --cds       mydata/cds_subset.fasta \
  --protein   mydata/proteins.fasta \
  --reference mydata/HK_reference.fasta \
  --outdir    custom_results
```

### Reuse existing results
You may configure the script to skip computation and only plot from existing CSVs (e.g., `results/tables/lovo_results.csv`).

### Debug mode (subset run)
Uncomment the sampling *if-block* around the debug section in the script (search for: `Debug: run only 12 human + 8 nonhuman species`).  
This samples **12 human** + **8 nonhuman** species for a quick run (~1.5 h). Re-comment the block for full runs.

---

## 7. Outputs

Under the chosen output directory (default `results/`), the pipeline writes:
- **tables/**: main label table, species summary, relabel stats, merged feature tables, `lovo_results.csv`, `train_test_results.csv`, feature-importance summaries.
- **figs/**: performance comparison, LOVO boxplots, UMAP projections, clustermaps, feature-importance plots.
- **logs/**: runtime logs.
- `pipeline_feature_set_coverage_summary.csv`.

---

## 8. Computational environment and reproducibility

All thresholds (e.g., `min_partial_length = 600 nt`), completeness filters, k-mer sizes, and class-weight options are set via script parameters.

Analyses were run in:
- Python ≥3.8
- Biopython ≥1.8
- pandas ≥1.5
- scikit-learn ≥1.2
- matplotlib ≥3.7
- seaborn ≥0.12
- umap-learn ≥0.5
- iFeatureOmegaCLI (current public release)

Randomised procedures used fixed seeds (42):
- Random Forest (`random_state=42`)
- Stratified 80/20 split (`random_state=42`)
- Group-aware 5-fold CV with shuffling (`StratifiedGroupKFold(..., random_state=42)`)

Intermediate outputs, logs, and provenance files are retained for reproducibility.  
On our hardware, a **full run (136 species)** with group-aware CV required ~13–14 h; the **debug subset (20 species)** ~1.5 h.

---

## 9. Citation

```bibtex
@software{jin_orthobunyavirus_pipeline_2025,
  author  = {Jin, Fengshi},
  title   = {Orthobunyavirus Host-Prediction Pipeline},
  year    = {2025},
  url     = {https://github.com/spike107/orthobunyavirus-host-prediction/releases/tag/v1.0.0},
  version = {v1.0.0}
}
```

---

## 10. License

MIT License (see `LICENSE`).  
Data are not distributed — users must obtain sequences and reference tables from their original sources (NCBI, ICTV, Virus-Host DB, etc.).

---

## 11. Contact

Please open GitHub issues or pull requests.
