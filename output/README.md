# 3-REAP Pipeline Output

**Project:** cmoore_020.bicli_pcf11  
**Pipeline:** 3-REAP — Paired-read Alternative Polyadenylation   
**Genome:** hg38  
**PAS reference:** PolyA_DB v4.1 (`hg38.PAS.main.tsv`)

---

## Samples and read counts

| Sample | Group | Trimmed R2 reads | Deduped R2 reads (BED) | PASS reads |
|---|---|---|---|---|
| Control_S1 | Control | 31,145,213 | 16,271,634 | 8,926,267 |
| Control_S2 | Control | 25,882,706 | 14,477,219 | 8,228,848 |
| Control_S3 | Control | 30,308,907 | 16,234,494 | 8,938,296 |
| Experimental_S4 | Experimental | 24,894,536 | 13,713,383 | 7,822,277 |
| Experimental_S5 | Experimental | 35,864,203 | 20,638,956 | 11,849,509 |
| Experimental_S6 | Experimental | 26,910,512 | 14,334,808 | 7,953,234 |

*Trimmed R2 reads*: reads surviving adapter and 5′-T-run trimming.  
*Deduped R2 reads*: read 2 alignments in the final BAM (R2 marks the 3′ cleavage end).  
*PASS reads*: LAPs (last aligned positions) assigned to a PolyA_DB v4.1 PAS within ±24 nt; used for all PASS bigWig tracks.

---

## Pipeline settings

| Parameter | Value |
|---|---|
| STAR version | 2.7.11b |
| Alignment: multi-map filter | 1 (uniquely mapping only) |
| Alignment: min score overlap | 0.2 |
| Soft-clip tolerance (`-misM`) | 2 nt |
| LAP–PAS assignment window | ±24 nt |
| Deduplication | Symlink to STAR BAM (UMI dedup disabled) |

---

## Directory contents

### `01_trim/` — Trimmed FASTQs

One subdirectory per sample. R1 and R2 were trimmed independently.

```
{sample}/
  {sample}_R1.fastq.gz                        R1: polyA, adapter, G-quality runs removed;
                                               reads with internal adapter (5′ orientation) discarded
  {sample}_R2_trimAdapt.5Ttrimmed.fastq.gz    R2: adapter trimmed, then 5′-T-run removed
                                               (trims homopolymer T-runs that indicate internal
                                               priming rather than true cleavage sites)
```

---

### `02_align/` — STAR alignment

One subdirectory per sample; paired-end reads were first repaired with BBMap `repair.sh` (~97–98% pair retention) before alignment.

```
{sample}/star/
  Aligned.sortedByCoord.out.bam     Coordinate-sorted BAM
  Aligned.sortedByCoord.out.bam.bai BAM index
  Log.final.out                     Per-sample alignment summary (unique %, multi-loci %, unmapped %)
  Log.out                           Full STAR run log
  Log.progress.out                  Real-time alignment progress log
  SJ.out.tab                        Splice junctions detected
  _STARtmp/BAMsort/                 Temporary sort files (can be deleted)

stats_STAR_summary.csv              Wide table: one row per sample, columns for total reads,
                                    uniquely mapped %, multi-loci %, and three unmapped categories
```

---

### `03_dedup/` — Deduplication

Currently a **symlink layer** — each `.dedup.bam` points directly to the STAR BAM in `02_align/`. UMI-based deduplication is disabled for this run (no UMI library prep).

```
{sample}/
  {sample}.dedup.bam    → symlink to 02_align/{sample}/star/Aligned.sortedByCoord.out.bam
```

---

### `04_bed/` — R2 BED files

Read 2 marks the RNA 3′ end (the cleavage and polyadenylation site). Only R2 alignments are carried forward from this step.

```
{sample}.dedup_R2.sorted.bed    R2 reads extracted (samtools -f 128 -F 4), converted to BED
                                 (bedtools bamtobed -cigar), sorted by chromosome
stats_deduped_reads.txt         Line counts per BED file (= deduped R2 read count per sample)
```

---

### `05_LAP_PAS/` — LAP definition and PAS assignment

LAP = **Last Aligned Position**: the 3′-most mapped base of each R2 alignment, used as a proxy for the cleavage site. LAPs within ±24 nt of a PolyA_DB v4.1 entry are called PASS (PolyA Site-Supporting) reads.

**Per-sample files:**

```
{sample}_LAP_bed.rds                R object (GRanges): all LAPs derived from the R2 BED
{sample}_LAP_polyAdb.rds            R object: LAPs with PolyA_DB distance annotations
{sample}_LAPs_24.PASS.rds           R object: PASS LAPs (≤24 nt from a PAS)
{sample}_stats.csv                  Per-sample LAP/PAS assignment statistics
{sample}_CIGAR_distrib.pdf          Soft-clipping distribution across reads
{sample}_LAP_polyAdb_distance.pdf   Distribution of LAP-to-PAS distances

{sample}_mapped_read.bed            All mapped R2 reads as BED (source for mapped_read bigWigs)
{sample}_PASS_bw_LAP24.bed          PASS reads within ±24 nt of any PAS (source for PASS_bw_LAP24 bigWigs)
{sample}_PASS_bw_LAP24_positon.bed  PASS reads collapsed to PAS position (PASS_bw_LAP24_positon bigWigs)
{sample}_PASS_bw_LAP24_PAS.bed      PASS reads at assigned PAS locus (PASS_bw_LAP24_PAS bigWigs)
```

**`PAS_quant/` — Quantification tables (primary downstream analysis input):**

```
{sample}_cluster.all.reads.csv    Per-sample PAS cluster read counts
cluster.all.reads.csv             Combined count matrix: rows = PAS clusters, columns = samples
                                  → primary input for differential APA analysis
stats_PAS.csv                     Combined per-sample assignment statistics across all samples
```

---

### `06_bigwig/` — Strand-specific bigWig tracks

CPM-normalized (counts per million), strand-specific coverage tracks ready for genome browser visualization (IGV, UCSC). Four track types, each in its own subdirectory.

**Track types:**

| Directory | Source BED | What it shows |
|---|---|---|
| `PASS_bw_LAP24/` | `*_PASS_bw_LAP24.bed` | PASS reads at LAPs within ±24 nt of a PAS |
| `PASS_bw_LAP24_positon/` | `*_PASS_bw_LAP24_positon.bed` | PASS reads collapsed to exact PAS position |
| `PASS_bw_LAP24_PAS/` | `*_PASS_bw_LAP24_PAS.bed` | PASS reads at assigned PAS locus |
| `mapped_read/` | `*_mapped_read.bed` | All uniquely mapped R2 reads (total 3′-end coverage) |

**File naming within each subdirectory:**

```
plus.{sample}.{type}.bw     Forward strand, CPM normalized
minus.{sample}.{type}.bw    Reverse strand, CPM normalized
mi.{sample}.{type}.bw       Reverse strand with negative signal values
                             (load mi.* for minus strand in IGV to display below the axis)
```

---

## Notes on `.dvc` files

Files with a `.dvc` extension are **DVC (Data Version Control) pointers**. They track the corresponding large file in remote storage and contain its checksum and size. The actual data file is the one without the `.dvc` suffix; the pointer file is used for version tracking and can be ignored during analysis.

---

## Key files for downstream analysis

| Purpose | File |
|---|---|
| Differential APA (primary counts) | `05_LAP_PAS/PAS_quant/cluster.all.reads.csv` |
| Per-sample PAS statistics | `05_LAP_PAS/PAS_quant/stats_PAS.csv` |
| Alignment QC | `02_align/stats_STAR_summary.csv` |
| Library sizes | `04_bed/stats_deduped_reads.txt` |
| Genome browser (all reads) | `06_bigwig/mapped_read/` |
| Genome browser (PASS reads) | `06_bigwig/PASS_bw_LAP24/` |
