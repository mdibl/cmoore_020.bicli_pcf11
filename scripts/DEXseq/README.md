# DEXseq.R — Differential Alternative Polyadenylation Analysis

**Script:** `scripts/DEXseq.R`
**Last updated:** 2026-06-25

---

## What this script does — and why

### The biological question

After the 3-REAP pipeline runs, we know how many RNA molecules from each sample ended at each poly-A site (PAS). But raw counts alone don't answer the biological question. What we want to know is: **does the relative usage of each PAS site change between the control and treatment conditions?**

Consider a gene that has two PAS: a proximal one (shorter 3' UTR) and a distal one (longer 3' UTR). If treatment causes cells to preferentially use the distal site, that gene is said to have undergone **3' UTR lengthening**. This shifts regulatory control — the longer 3' UTR exposes binding sites for microRNAs and RNA-binding proteins that weren't accessible before, potentially changing how much protein is made, how stable the mRNA is, or where in the cell it goes.

The fraction of reads from a gene that end at a specific PAS is called **PSI (Percent Site Inclusion)** or **fractional usage**. A PSI of 0.7 at the distal site means 70% of that gene's poly-A reads end at the distal PAS. We test whether PSI changes significantly between conditions.

### The statistical method: DEXSeq

This script uses **DEXSeq**, a Bioconductor package originally developed to detect differential exon usage in RNA-seq data. We repurpose it for APA by treating each PAS as if it were an exon within its gene. DEXSeq's core insight is that comparing *relative* usage (each PAS as a fraction of its gene's total counts) is statistically more powerful than comparing absolute counts, because it controls for variation in overall gene expression.

The DEXSeq model tests each PAS independently with the design `~ sample + exon + condition:exon`:
- `sample` absorbs per-sample baseline differences (the blocking factor)
- `condition:exon` is the actual test: does the Treatment condition change the fractional usage of this PAS?

A significant result (adjusted p-value < threshold) means the relative usage of that PAS is statistically different between control and treatment, after accounting for how many reads the gene received overall.

### APA direction: Lengthened vs. Shortened

For each significant PAS, the script classifies the direction of the change:

- **Distal PAS**: the PAS with the most extreme genomic position in the 3' direction (highest coordinate on + strand; lowest coordinate on − strand). This corresponds to the longest possible 3' UTR.
- **Lengthened**: the distal PAS gains relative usage in treatment → longer 3' UTRs
- **Shortened**: a proximal PAS gains relative usage in treatment → shorter 3' UTRs
- **Mixed**: a gene has some PAS going up and some going down

---

## Quick-start: how to run the script

### 1. Edit the configuration block

Open `DEXseq.R`. Change only the `USER CONFIGURATION — EDIT THESE` section at the top:

```r
BASE_DIR    <- "/path/to/your/project"          # project root
COUNTS_CSV  <- "path/to/cluster.all.reads.csv"          # PAS count matrix from 3-REAP
SAMPLES_TXT <- "path/to/sample_name.txt"                # sample list (or set to NULL)
ANNO_FILE   <- "data/PolyA_DB_v4.1/hg38.PAS.main.tsv"  # same file used by 3-REAP
DESIGN_FILE <- "path/to/design.txt"                     # experimental design table

CTRL_LABEL  <- "Control"     # must match values in design file's 'condition' column
TRTMT_LABEL <- "Treatment"   # must match values in design file's 'condition' column

GROUPING_VAR <- NULL   # NULL = run one analysis on all samples
                       # set to "timepoint" (or any colData column) to run per level

GENES_OF_INTEREST <- c("GENE1", "GENE2")  # genes to plot PSI bar charts for

USE_RUV <- FALSE   # set to TRUE to enable RUVseq batch correction (see below)
RUV_K   <- 1       # number of unwanted variation factors to estimate (1 is usually enough)
```

### 2. Prepare input files

**Count matrix** (`cluster.all.reads.csv`): produced by 3-REAP step 6. A CSV where the first column is `hit_PAS_ID` and the remaining columns are samples.

**PAS annotation file** (`data/PolyA_DB_v4.1/hg38.PAS.main.tsv`): use the **same file** that was provided to the 3-REAP pipeline. The `hit_PAS_ID` values in the count matrix are taken directly from this database, so the IDs will match. The script auto-detects the PolyA_DB format:
- v3.2: explicit `Chromosome`, `Position`, `Strand` columns
- v4.1: coordinates encoded in `PAS_ID` as `chr:strand:position` (auto-parsed; default for this project)

**Design file** (`design.txt`): a tab-delimited file, first column = sample names (matching column names in the count matrix), remaining columns = metadata. Must include at minimum:

```
sample_name    condition
Control_rep1   Control
Control_rep2   Control
Treatment_rep1 Treatment
Treatment_rep2 Treatment
```

If you have a timepoint or batch column, include it here. Set `GROUPING_VAR <- "timepoint"` to run a separate DEXSeq analysis per timepoint level.

**Sample name file** (`sample_name.txt`, optional): a plain text file listing sample names, one per line. If set to NULL, all non-`hit_PAS_ID` columns in the count matrix are used.

### 3. Run

```r
Rscript scripts/DEXseq.R
```

Or open in RStudio and run interactively. The script prints progress messages for each step.

---

## Software requirements (R packages)

Install via Bioconductor/CRAN if missing:

```r
BiocManager::install(c("DEXSeq", "DESeq2", "GenomicRanges", "SummarizedExperiment"))
install.packages(c("ggplot2", "ggrepel", "scales", "limma", "dplyr",
                   "tidyr", "tibble", "matrixStats", "stringr"))
```

---

## Understanding the analysis: step by step

### Filtering
Before testing, PAS sites are filtered to retain only:
1. **3' UTR PAS** — sites in the 3' untranslated region (set by `PAS_TYPE_REGEX`). Intronic and other PAS are excluded because APA in 3' UTRs is the biologically relevant type for most questions.
2. **Sufficient read coverage** — PAS with fewer than `MIN_TOTAL_READS` supporting reads in a group are removed. This prevents testing sites where we have too little evidence.
3. **Genes with ≥2 PAS** — DEXSeq tests relative usage within genes; a gene with only one observable PAS cannot have differential isoform usage.

### Normalization
DEXSeq estimates a **size factor** for each sample — a scaling factor that accounts for differences in sequencing depth (some samples have more reads than others). All counts are internally divided by their sample's size factor before testing.

### Testing
DEXSeq fits a statistical model to each PAS and tests whether the Treatment:exon interaction term is significantly non-zero — i.e., whether the fractional usage of that PAS changes in treatment vs. control. This produces a p-value per PAS. Multiple testing correction (Benjamini-Hochberg FDR) is applied to give an adjusted p-value (`padj`).

Additionally, the Lancaster procedure (`perGeneQValue`) combines per-PAS p-values into a per-gene q-value that summarizes whether *any* PAS in the gene shows significant differential usage.

### Fold changes
DEXSeq reports a **log2 fold change** for each PAS: the log2 ratio of normalized fractional usage (Treatment / Control). However, for biological interpretation we report **ΔPSI** (Treatment PSI − Control PSI) in the plots, because it is directly interpretable as a percentage-point change in isoform usage.

A ΔPSI of +0.20 at the distal PAS means the distal site gained 20 percentage points of relative usage in treatment.

---

## Output files

All outputs are written under `OUT_BASE` (set in configuration). Structure:

```
{OUT_BASE}/
├── pas_results[.group].csv     ← per-PAS DEXSeq results with APA direction
├── pas_usage[.group].csv       ← per-PAS mean PSI + SEM per condition
├── gene_summary[.group].csv    ← per-gene summary
└── plots/
    ├── library_sizes.png
    ├── PCA_global.png                  ← (USE_RUV=FALSE only)
    ├── PCA_global_uncorrected.png      ← (USE_RUV=TRUE only)
    ├── PCA_global_RUV_corrected.png    ← (USE_RUV=TRUE only; visual only)
    ├── RUV_W_factors.png               ← (USE_RUV=TRUE only)
    ├── PCA_normalized[.group].png
    ├── size_factors[.group].png
    ├── dispersion[.group].png
    ├── pvalue_hist[.group].png
    ├── MA[.group].png
    ├── volcano[.group].png
    └── usage/
        └── {GENE}[.group].png   ← PSI bar charts, one file per gene
```

`[.group]` is appended only when `GROUPING_VAR` is set (e.g., `.3d` for a "3day" timepoint level). With `GROUPING_VAR <- NULL`, files have no suffix: `pas_results.csv`, `gene_summary.csv`, etc.

---

### `pas_results.csv` — per-PAS results

The main results table. One row per PAS site that passed filtering and was tested.

| Column | Meaning |
|---|---|
| `groupID` | Gene identifier (Gene Symbol if available, else Ensembl/RefSeq) |
| `featureID` | PAS site ID (from PolyA_DB) |
| `gene` | Same as groupID |
| `PAS_type` | Annotation category (e.g., `3'UTR`, `intron`) |
| `APA_direction` | `Lengthened`, `Shortened`, `Ambiguous`, or `NA` |
| `exonBaseMean` | Mean normalized read count at this PAS across all samples |
| `pvalue` | DEXSeq per-PAS p-value |
| `padj` | FDR-adjusted p-value (Benjamini-Hochberg) |
| `log2fold_Treatment_v_Control` | log2(Treatment fractional usage / Control fractional usage) |
| `countData.Control_rep1`, ... | Raw read counts per sample |
| `genomicData.seqnames/start/strand` | Genomic coordinates of the PAS |

**How to identify significant APA events:** filter for `padj < 0.05` (or your threshold). Further filter by `APA_direction == "Lengthened"` or `"Shortened"` to focus on directional 3' UTR changes.

---

### `gene_summary.csv` — per-gene summary

One row per gene. Summarizes across all tested PAS within that gene.

| Column | Meaning |
|---|---|
| `groupID` | Gene identifier |
| `n_PAS` | Total PAS tested for this gene |
| `n_sig` | Number of PAS with `padj < PADJ_CUT` |
| `n_big` | Number of significant PAS with `\|log2FC\| >= LFC_CUT` |
| `min_padj` | Most significant adjusted p-value among this gene's PAS |
| `max_absL2FC` | Largest absolute log2 fold change among significant PAS |
| `dir_consensus` | `Lengthened_only`, `Shortened_only`, or `Mixed` |
| `perGeneQ` | Lancaster-combined gene-level q-value (primary ranking metric) |

**How to read this table:** sort by `perGeneQ` ascending. Genes at the top have the strongest overall evidence for differential APA. `dir_consensus` tells you whether the change is directional (one dominant direction) or complex (multiple PAS moving in different directions, which can indicate alternative internal promoter usage or complex regulatory changes).

---

### `pas_usage.csv` — PSI values

Same rows as the per-PAS results file, plus four additional columns:

| Column | Meaning |
|---|---|
| `meanUsage_Control` | Mean PSI across control replicates (0 to 1; 0.5 = 50% of gene's reads) |
| `meanUsage_Treatment` | Mean PSI across treatment replicates |
| `seUsage_Control` | Standard error of PSI across control replicates |
| `seUsage_Treatment` | Standard error of PSI across treatment replicates |

PSI (fractional usage) = reads at this PAS / total reads at all PAS of this gene, per sample, then averaged. This is the most biologically interpretable metric: a PAS with `meanUsage_Control = 0.25` accounts for 25% of that gene's poly-A reads in controls.

---

## Quality control plots

All QC plots are saved to `plots/`.

### `library_sizes.png`
Bar chart of total PAS read counts per sample (before any filtering). Bars should be roughly similar height across samples. Large differences (>3×) may indicate a failed sample or a sample that needs closer attention during normalization. This is assessed before the analysis runs.

### `PCA_global.png` / `PCA_global_uncorrected.png`
Principal component analysis on all samples simultaneously, using variance-stabilized 3' UTR PAS counts. Points that cluster by condition (e.g., all Controls together) indicate that the experimental effect is larger than technical variation — a good sign. Points that cluster by batch or unexpected metadata suggest confounding. When `USE_RUV=FALSE` this file is named `PCA_global.png`; when `USE_RUV=TRUE` it is named `PCA_global_uncorrected.png` so the pre-correction structure is preserved for comparison.

### `PCA_global_RUV_corrected.png` *(USE_RUV=TRUE only)*
The same PCA after regressing the estimated W factors out of the VST matrix using `limma::removeBatchEffect`. This is **for visualization only** — it shows what the data look like after removing the batch signal, but does not affect the statistical model. Compare this to the uncorrected PCA to see how much the W factors explained.

### `RUV_W_factors.png` *(USE_RUV=TRUE only)*
Bar chart of the estimated RUV W factor(s) per sample, colored by condition. Samples that drove the batch factor will have large positive or negative bars. Ideally this factor should not be correlated with condition — if it is, the correction may over-remove true signal.

### `PCA_normalized[.group].png`
PCA on the post-normalization counts within each analysis group. Similar interpretation as above, but after DEXSeq's size factor normalization.

### `size_factors[.group].png`
Bar chart of the normalization factors assigned by DEXSeq to each sample. All bars should be close to 1.0 (the dashed reference line). A sample with a very high or low size factor (>2× away from 1) may be an outlier.

### `dispersion[.group].png`
This plot shows the estimated variance of each PAS's counts across samples. DEXSeq fits a curve (shown in red) through the cloud of per-PAS dispersion estimates (black dots) and uses that curve to smooth individual estimates toward the genome-wide trend — a technique called shrinkage. Well-behaved data has a tight cloud around the fitted curve with a mean that decreases at higher mean counts.

### `pvalue_hist[.group].png`
Histogram of raw (unadjusted) p-values. The dashed red line shows what a uniform distribution would look like. A well-powered experiment with true effects shows a spike of small p-values near 0 (real signal) and a roughly flat distribution from 0.1–1.0 (null PAS). A completely flat or anti-conservative (U-shaped) histogram suggests model misspecification.

---

## Diagnostic plots

Saved to `plots/`.

### `MA_plot[.group].png`
Plots log2 fold change (Y) versus mean count (X) for every tested PAS. Significant hits (padj < threshold) are shown in red. A symmetric cloud centered at 0 with random scatter indicates no systematic bias. Up/down counts are annotated in the corners.

### `volcano[.group].png`
Plots log2 fold change (X) versus −log10(padj) (Y). Points in the upper corners are the most significant and most changed PAS. The top hits by significance are labeled with gene names. Color categories:
- **Red** — significant AND large fold change (`padj < threshold` AND `|log2FC| >= LFC_CUT`)
- **Gold** — significant but small fold change
- **Gray** — not significant

---

## PSI bar charts (relative usage plots)

Saved to `plots/usage/`, one PNG per gene per group.

One plot per gene in `GENES_OF_INTEREST`, per analysis group. Each plot shows:

- **X axis:** PAS sites ordered from 5' to 3' along the transcript (proximal on left, distal on right). The distal PAS is labeled `(distal)`.
- **Y axis:** Mean fractional usage (PSI), expressed as a percentage. The values for all PAS in the gene add up to 100%.
- **Bars:** Side-by-side bars for Control (left) and Treatment (right) at each PAS.
- **Error bars:** ±1 standard error of the mean (SEM) across replicates within each condition. Larger error bars mean more variability between replicates.
- **Labels:** Each PAS is labeled with **ΔPSI** (Treatment PSI − Control PSI, as a decimal fraction). A positive ΔPSI means more usage in treatment. Significance markers are appended: `*` (padj < 0.05), `**` (padj < 0.01), `***` (padj < 0.001).

**How to read these plots for a collaborator:**

Look at whether the bars shift toward the right (distal) or left (proximal) in the treatment vs. control. If the rightmost bar (distal PAS, labeled `(distal)`) gets taller in treatment, the gene's 3' UTR is getting *longer* in treatment. If a left bar (proximal PAS) gets taller while the distal bar shrinks, the 3' UTR is getting *shorter*. The ΔPSI value tells you the magnitude of that shift in percentage points.

---

## Interpreting the results: a guide for collaborators

**What does it mean when a gene shows differential APA?**

It means cells under the two conditions preferentially choose different poly-A sites for that gene. This changes the length of the mRNA's 3' UTR without necessarily changing how much of the protein is made (though it can do that too).

A longer 3' UTR exposes more regulatory sequences — microRNA binding sites, RNA-binding protein sites, adenylation signals — that can reduce mRNA stability or protein output. A shorter 3' UTR "hides" those sites. APA is therefore a mechanism for fine-tuning gene regulation without changing the gene's transcription rate.

**Which genes should I pay attention to?**

The gene-summary table (`gene_summary.csv`) ranked by `perGeneQ` is the best place to start. Genes with:
- Low `perGeneQ` (close to 0) — strong statistical evidence for differential APA
- `dir_consensus = "Lengthened_only"` or `"Shortened_only"` — clean, directional changes
- High `n_sig` — multiple PAS in the gene are changing, not just one

**What does a significant padj mean vs. a large ΔPSI?**

`padj` tells you how confident we are that the change is real (statistical significance). ΔPSI tells you how large the change is (biological effect size). A small ΔPSI can be highly significant if we have many replicates. A large ΔPSI may not reach significance if the data are noisy or replicates are few. Both matter: ideally you want both low padj and large |ΔPSI|.

**What is "perGeneQ"?**

The per-gene q-value (Lancaster method) combines the evidence from all PAS sites within a gene into a single gene-level significance score. It is more conservative than just taking the minimum padj, but it accounts for the fact that testing multiple PAS per gene inflates the chance of at least one false positive. Use `perGeneQ` for gene-level prioritization.

---

## Running for different experimental designs

### Simple two-condition (no grouping)
Leave `GROUPING_VAR <- NULL`. All samples in the design file are analyzed together.

### Multiple timepoints (or batches as separate analyses)
Set `GROUPING_VAR <- "timepoint"` (or whatever the column is called in your design file). The script automatically detects all levels in that column and runs a separate DEXSeq analysis for each — for example, one for day 3 samples only, and one for day 5 samples only. Output files are suffixed accordingly (`.3d`, `.5d`).

### Adding a batch covariate (blocking factor)
DEXSeq's `~ sample + exon + condition:exon` model absorbs all per-sample gene-level variation through the `sample` term, which handles differences in overall expression between samples. However, it does **not** correct for batch effects on *relative* PAS usage (PSI). If replicates cluster by batch rather than condition in the PCA, use RUVseq.

### RUVseq batch correction (`USE_RUV <- TRUE`)
Set `USE_RUV <- TRUE` when your PCA shows replicates clustering by something other than condition — for example, if one replicate from each condition groups separately from the others, suggesting a technical batch effect.

RUVseq uses the within-condition replicate structure to estimate "unwanted variation" W factors — sources of variability that are shared within replicates of the same condition and therefore cannot be biological. `RUV_K` controls how many such factors to estimate (1 is almost always sufficient; try 2 if the first factor does not resolve the outlier clustering).

When `USE_RUV=TRUE` the DEXSeq design becomes:
```
~ sample + exon + W_1:exon + condition:exon
```
The `W_1:exon` term tests whether the batch factor changes relative PAS usage fractions within genes — this is the dimension that `sample` alone cannot capture. The `condition:exon` term retains the same interpretation as before.

**Important:** RUVseq corrects for unwanted variation that exists consistently across replicates. It cannot rescue experiments where replicates are fundamentally incomparable. Always compare the uncorrected and corrected PCA plots before trusting the RUV-corrected results.

---

## Key thresholds and what they control

| Parameter | Default | Effect |
|---|---|---|
| `PAS_TYPE_REGEX` | `^3'UTR` | Which PAS types to include. Change to `""` to include all types. |
| `MIN_TOTAL_READS` | 10 | Minimum total reads at a PAS in the group. Lower = more PAS tested, more noise. |
| `MIN_PER_CONDITION` | 1 | Minimum samples with >0 reads per condition. Ensures both conditions have detectable signal. |
| `PADJ_CUT` | 0.05 | FDR threshold. Lower = fewer but more confident hits. |
| `LFC_CUT` | 0.5 | log2 fold change threshold for "big" hits in the gene summary and volcano plot coloring. |
| `NTOP_PCA` | 2000 | Number of most-variable PAS to use for PCA. 500–5000 is typical. |

---

## Notes on the method

**Why DEXSeq and not a simpler test?**

DEXSeq's key advantage is that it explicitly models the within-gene relative usage using a negative binomial model, which handles the overdispersion typical of count data. The `sample` term absorbs all inter-sample variability, making the test highly specific to isoform usage changes rather than overall expression changes. Simpler approaches (e.g., Fisher's exact test on counts) do not account for between-replicate variance and over-report false positives.

**Why PSI and not log2FC for the bar plots?**

The DEXSeq log2 fold change is a ratio-of-ratios on the count scale (internally: Treatment PSI / Control PSI, log2-transformed). While it is useful for volcano/MA plots, it is not intuitively interpretable. ΔPSI (simply Treatment PSI − Control PSI) is expressed in the same units as the bars themselves and is directly relatable to the plot — a ΔPSI of +0.15 means the bar got 15 percentage points taller.
