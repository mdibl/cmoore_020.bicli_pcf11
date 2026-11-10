# DEXseq_allreps — All 6 Replicates, Standard Run

**Run date:** 2026-06-25
**Script:** `scripts/DEXseq/DEXseq.R`
**Output base:** `output/DEXseq_allreps/`

---

## What this run was

A standard DEXSeq APA analysis with no batch correction and no sample exclusions. All 6 samples were included as-is.

**Samples:**
- Control rep1, rep2, rep3
- Experimental rep1, rep2, rep3

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

## Purpose

Baseline result. All replicates included; no adjustments made. Use this run to compare against `DEXseq_RUVs` (same samples, with batch correction) to see how much the batch factor affected the results.

Note: PCA showed rep3 in both conditions clustering separately from rep1/rep2, suggesting a replicate-level batch effect. This run does **not** correct for that. See `DEXseq_RUVs` for the corrected version.

---

## Output files

```
output/DEXseq_allreps/
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
