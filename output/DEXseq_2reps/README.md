# DEXseq_2reps — 2 Replicates Per Condition (rep3 excluded)

**Run date:** 2026-06-25
**Script:** `scripts/DEXseq/DEXseq.R`
**Output base:** `output/DEXseq_2reps/`

---

## What this run was

DEXSeq APA analysis restricted to rep1 and rep2 from each condition. Rep3 was excluded from both Control and Experimental groups by editing `data/sample_name.txt` and `data/design.txt` to list only the 4 retained samples.

**Samples:**
- Control rep1, Control rep2
- Experimental rep1, Experimental rep2

**Key settings (all others at default):**

| Parameter | Value |
|---|---|
| `USE_RUV` | `FALSE` — no batch correction |
| `CTRL_LABEL` | `Control` |
| `TRTMT_LABEL` | `Experimental` |
| `PAS_TYPE_REGEX` | `^3'UTR` |
| `MIN_TOTAL_READS` | 10 |
| `PADJ_CUT` | 0.05 |

**DEXSeq design:** `~ sample + exon + condition:exon`

---

## Why rep3 was excluded

PCA on the full 6-sample dataset showed rep3 from both conditions clustering separately from rep1/rep2 — a pattern consistent with a technical batch effect introduced during sample preparation or sequencing of those two samples. Removing both rep3s eliminates that variation entirely at the cost of reduced statistical power (n=2 per condition instead of n=3).

---

## Purpose

Sensitivity check. Results from this run can be compared to `DEXseq_RUVs` (which retains rep3 but corrects for the batch effect) to assess whether:
- Hits found in both runs are robust — they survive both approaches
- Hits unique to `DEXseq_RUVs` are likely real signal recovered by the correction
- Hits unique to this run may reflect false positives from the reduced sample size

With only 2 replicates per condition, statistical power is limited — expect fewer significant hits and wider confidence intervals than the 3-replicate runs.

---

## Output files

```
output/DEXseq_2reps/
├── pas_results.csv       ← per-PAS DEXSeq results + APA direction
├── pas_usage.csv         ← per-PAS mean PSI + SEM per condition
├── gene_summary.csv      ← per-gene summary ranked by perGeneQ
└── plots/
    ├── library_sizes.png
    ├── PCA_global.png
    ├── PCA_normalized.png
    ├── size_factors.png
    ├── dispersion.png
    ├── pvalue_hist.png
    ├── MA.png
    ├── volcano.png
    └── usage/
        └── PCF11.png, TAB2.png, ICAM1.png
```

For column definitions see `scripts/DEXseq/README.md`.
