# NBa-SEic Pipeline 🧬  
*Generate non-binary super-enhancer tables from ROSE peak files*

---

## 1 Quick start

```bash
# 1) Create and activate the conda environment
conda env create -f environment.yml
conda activate nbase-seic

# 2) Run the pipeline on the bundled tiny test-set
./run_pipeline.sh example_data results
```

Four CSV tables will appear in `results/Tables/`:

| File                | Brief contents                                                         |
|---------------------|------------------------------------------------------------------------|
| **Set2__Table1.csv** | Consolidated SE loci (12.5 kb merge) + binary presence flags          |
| **Set2__Table2.csv** | All SE/TE elements per locus + overlap weights                        |
| **Set2__Table3.csv** | Locus activity per sample (peak & weighted) + SE/TE counts            |
| **Set2__Table4.csv** | Feature matrix (median-norm → log1p → z-score) + binary flags        |

---

### Run on your own data

```bash
./run_pipeline.sh /path/to/rose_peaks  my_results
```

`/path/to/rose_peaks` must contain **two** ROSE-generated BED files for **each** sample:

| File pattern                  | Description                              |
|-------------------------------|------------------------------------------|
| `SE_<sample>_SE_mm10.bed`     | Super-enhancer peaks (case vs. input)    |
| `SE_<sample>_TE_mm10.bed`     | Typical-enhancer peaks                   |

The script is species-agnostic—only the filenames need to match these patterns.

---

## 2 What the pipeline does

1. **Read ROSE peaks**  
   Load SE & TE BED files for each sample.  

2. **Consolidate loci**  
   Merge overlapping SEs within 12 500 bp to define shared loci.  

3. **Map enhancers**  
   Intersect consolidated loci with the union of SE + TE elements.  

4. **Compute activity**  
   Calculate peak and weighted H3K27ac signal per locus & sample.  

5. **Engineer features**  
   Median-normalization → global imputation (1 × 10⁻⁶) → log1p → z-score.  

6. **Export tables**  
   Write the four ready-to-analyse CSV files.

> **Note:** This pipeline runs *on top of* ROSE—it does *not* process raw FASTQ files, so you maintain full control of your peak-calling parameters.

---

## 3 Requirements

| Software                                     | Installed via                                |
|----------------------------------------------|----------------------------------------------|
| Python 3.11, pandas, numpy, scikit-learn      | `conda env create -f environment.yml`         |
| bedtools ≥ 2.30                              | Included in the same conda environment       |
| Linux/macOS shell *or* Windows + Git Bash/WSL | Required for `sort` & `bedtools` commands     |

---

## 4 Troubleshooting

| Symptom                                         | Fix                                                                 |
|-------------------------------------------------|---------------------------------------------------------------------|
| `bedtools: command not found`                   | Activate the conda environment: `conda activate nbase-seic`.        |
| *No such file* `SE_<sample>_SE_mm10.bed`        | Verify your filenames and the `--data_dir` path.                    |
| Memory errors on very large cohorts             | Split samples by tissue or chromosome and run in smaller batches.   |

---

## 5 Citation

If you use this pipeline, please cite:

> **Osintseva E.D. _et al._**  
> _A Non-Binary Approach to Super-Enhancer Identification and Clustering: A Dataset for Tumor and Treatment-Associated Dynamics in Mouse Tissues._  
> **Data** (2025).
