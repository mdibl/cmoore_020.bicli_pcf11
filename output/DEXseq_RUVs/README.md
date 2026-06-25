# DEXseq_RUVs — All 6 Replicates with RUVseq Batch Correction

**Run date:** 2026-06-25
**Script:** `scripts/DEXseq/DEXseq.R`
**Output base:** `output/DEXseq_RUVs/`

---

## What this run was

DEXSeq APA analysis on all 6 samples with RUVseq batch correction enabled (`USE_RUV=TRUE`, `RUV_K=2`). This run addresses the rep3 outlier effect visible in PCA, where rep3 from both conditions clustered away from rep1/rep2.

**Samples:**
- Control rep1, rep2, rep3
- Experimental rep1, rep2, rep3

**Key settings:**

| Parameter | Value |
|---|---|
| `USE_RUV` | `TRUE` |
| `RUV_K` | `2` — two unwanted variation factors estimated |
| `CTRL_LABEL` | `Control` |
| `TRTMT_LABEL` | `Experimental` |
| `PAS_TYPE_REGEX` | `^3'UTR` |
| `MIN_TOTAL_READS` | 10 |
| `PADJ_CUT` | 0.05 |

**DEXSeq design:** `~ sample + exon + W_1:exon + W_2:exon + condition:exon`

The `W_1:exon` and `W_2:exon` terms correct for batch effects on *relative* PAS usage (PSI) within genes. They are estimated by RUVseq using variation among within-condition replicates — any variation shared between rep3 Control and rep3 Experimental but not present in rep1/rep2 is captured in W and removed from the test.

---

## Purpose

The preferred run for this dataset. Retains all 3 replicates per condition (greater statistical power than `DEXseq_2reps`) while controlling for the technical variation causing rep3 to cluster separately. Compare results to `DEXseq_allreps` to see the effect of the correction, and to `DEXseq_2reps` to confirm that the correction is recovering real signal rather than masking noise.

---

## Output files

```
output/DEXseq_RUVs/
├── pas_results.csv       ← per-PAS DEXSeq results + APA direction
├── pas_usage.csv         ← per-PAS mean PSI + SEM per condition
├── gene_summary.csv      ← per-gene summary ranked by perGeneQ
└── plots/
    ├── library_sizes.png
    ├── PCA_global_uncorrected.png   ← raw PCA showing batch structure
    ├── PCA_global_RUV_corrected.png ← visual-only; W factors regressed out
    ├── RUV_W_factors.png            ← per-sample W_1 and W_2 values
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
For RUVseq interpretation see the "RUVseq batch correction" section of that README.
