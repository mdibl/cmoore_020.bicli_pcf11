### 3-REAP-paired reads


#####################################################
### USER-DEFINED CONFIGURATION — EDIT THESE! ###
#####################################################

script_path=/compbio/scratch/cnobrega/3_REAP/scripts
fastq_path=/compbio/data/ClaireMoore/cmoore_019.admera
work_path=/compbio/scratch/cnobrega/3_REAP/cmoore_019.admera/output
sample_file=/compbio/scratch/cnobrega/3_REAP/cmoore_019.admera/sample.txt

threads_per_job=2 # cores for STAR
max_jobs=4 # number of samples to process in parallel

star_index_path=/compbio/scratch/cnobrega/3_REAP/ref/star_index_hg38 # for Step 3, STAR index. # If you have no star index, please build one first.
refPAS_file=/compbio/scratch/cnobrega/3_REAP/ref/human.PAS.hg38.txt # for Step 6, PAS annotation file from PolyA_DBv3.2
chromsizes=/compbio/scratch/cnobrega/3_REAP/ref/hg38.chrom.sizes # for Step 7, Chromosome sizes file

echo "Starting 3-REAP pipeline..."
#=============================================================
### step 1. UMI extraction
#=============================================================

#echo "Step 1: UMI extraction"
#
#cd $fastq_path
#while IFS=$'\t' read -r sample fastq_1 fastq_2; do
#	outpath=${work_path}/s1_fastq_umi/${sample}
#	mkdir -p ${outpath}
#	umi_tools extract --extract-method=string --bc-pattern X --bc-pattern2 NNNNNNNNNN \
#	-I $fastq_1 \
#	--read2-in=$fastq_2 \
#	-L ${outpath}/extraction_log \
#	-S ${outpath}/${sample}_R1.fastq.gz \
#	--read2-out=${outpath}/${sample}_R2.fastq.gz
#done < $sample_file

#=============================================================
### step 2. Trimming adapters and polyA/polyT sequences
#=============================================================
echo "Step 2: Adapter and polyA/polyT trimming"

cd $fastq_path

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
	{
		outpath=${work_path}/s2_fastq_cutadapt/${sample}
		outfile_r1=${outpath}/${sample}_R1.fastq.gz
		if [[ -f "$outfile_r1" ]]; then
			echo "Skipping $sample read 1 (already trimmed)"
			exit 0
		fi
		
		echo "Processing $sample read 1..."
		mkdir -p ${outpath}
		cutadapt --quiet -j 24 -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 $fastq_1 | \
		cutadapt --quiet -j 24 -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
		cutadapt --quiet -j 24 -m 20 -O 3 -a "r1polyA=A{18}" - | \
		cutadapt --quiet -j 24 -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" -o ${outpath}/${sample}_R1.fastq.gz --discard-trimmed -
	} & 
	((job_count++))
	if ((job_count >= max_jobs)); then
    	wait -n
    	((job_count--))
    fi
done < $sample_file
wait

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    {
		outpath=${work_path}/s2_fastq_cutadapt/${sample}
		outfile_r2=${outpath}/${sample}_R2_trimAdapt.5Ttrimmed.fastq.gz
		if [[ -f "$outfile_r2" ]]; then
			echo "Skipping $sample read 2 (already trimmed)"
			exit 0
		fi
		
		echo "Processing $sample read 2..."
		outpath=${work_path}/s2_fastq_cutadapt/${sample}
		mkdir -p ${outpath}
		cutadapt --quiet -j 24 -m 18 -O 3 --nextseq-trim=10 -a "rd2Adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=3;max_error_rate=0.100000" $fastq_2 | \
		cutadapt --quiet -j 24 -m 18 -O 20 -g "rd2Adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=20" -o ${outpath}/${sample}_R2_trimAdapt.fastq.gz --discard-trimmed -
		### Reverse reads contain 5'Ts, remove remaining 5'Ts before the alignment
		python ${script_path}/trim_N-5T_fq.py --rawfastq ${outpath}/${sample}_R2_trimAdapt.fastq.gz --out_dir ${outpath} --random_NT_len 0
    } & 
    ((job_count++))
	if ((job_count >= max_jobs)); then
    	wait -n
    	((job_count--))
    fi
done < $sample_file
wait

#=============================================================
### step 3. Genome alignment using paired reads
#=============================================================
# Since the fastq files were preprocessed separately, they need to be repaired to re-establish read pairing before alignment. repair.sh is a tool of bbmap.
# If you have no star index, please build one first.
echo "Step 3: Repair read pairs and align with STAR"

mkdir -p ${work_path}/s2_fastq_cutadapt_pair

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
  {
    repaired_r1=${work_path}/s2_fastq_cutadapt_pair/${sample}/${sample}_R1.fastq.gz
    repaired_r2=${work_path}/s2_fastq_cutadapt_pair/${sample}/${sample}_R2.fastq.gz
    star_bam=${work_path}/s3_star_align/${sample}/Aligned.sortedByCoord.out.bam

    if [[ -f "$repaired_r1" && -f "$repaired_r2" && -f "$star_bam" ]]; then
      echo "Skipping $sample (already repaired and aligned)"
      exit 0
    fi
    echo "Starting: $sample"
    # Create necessary output directories
    mkdir -p ${work_path}/s2_fastq_cutadapt_pair/${sample}
    mkdir -p ${work_path}/s3_star_align/${sample}

    # Repair
    repair.sh overwrite=t \
      in1=${work_path}/s2_fastq_cutadapt/${sample}/${sample}_R1.fastq.gz \
      in2=${work_path}/s2_fastq_cutadapt/${sample}/${sample}_R2_trimAdapt.5Ttrimmed.fastq.gz \
      out1=${work_path}/s2_fastq_cutadapt_pair/${sample}/${sample}_R1.fastq.gz \
      out2=${work_path}/s2_fastq_cutadapt_pair/${sample}/${sample}_R2.fastq.gz \
      outs=${work_path}/s2_fastq_cutadapt_pair/${sample}/${sample}_unpaired.fastq.gz \
      repair

    # Align
    STAR --runThreadN $threads_per_job --readFilesCommand zcat --genomeDir ${star_index_path} \
      --readFilesIn ${work_path}/s2_fastq_cutadapt_pair/${sample}/${sample}_R1.fastq.gz ${work_path}/s2_fastq_cutadapt_pair/${sample}/${sample}_R2.fastq.gz \
      --outReadsUnmapped None \
      --outFilterType BySJout --limitBAMsortRAM 4000000000 \
      --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${work_path}/s3_star_align/${sample}/ \
      --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.2 --outFilterMatchNminOverLread 0.2

  } &
  ((job_count++))
  if ((job_count >= max_jobs)); then
    wait -n
    ((job_count--))
  fi
done < $sample_file
wait

## Statistical alignment rate
echo "Compiling STAR alignment statistics..."
cd ${work_path}/s3_star_align
output_file="stats_STAR_summary.csv"
echo "Sample,Total Reads,Uniquely Mapped Reads (%),Mapped to Multiple Loci (%),Unmapped: Too Many Mismatches (%),Unmapped: Too Short (%),Unmapped: Other (%)" > $output_file
for logfile in */*Log.final.out; do
    #sample=$(basename $(dirname $logfile))
    sample=${logfile%%Log.final.out}
    total_reads=$(grep "Number of input reads" $logfile | awk '{print $NF}')
    uniquely_mapped=$(grep "Uniquely mapped reads %" $logfile | awk '{print $NF}')
    multiple_mapped=$(grep "% of reads mapped to multiple loci" $logfile | awk '{print $NF}')
    unmapped_mismatches=$(grep "% of reads unmapped: too many mismatches" $logfile | awk '{print $NF}')
    unmapped_short=$(grep "% of reads unmapped: too short" $logfile | awk '{print $NF}')
    unmapped_other=$(grep "% of reads unmapped: other" $logfile | awk '{print $NF}')
    echo "$sample,$total_reads,$uniquely_mapped,$multiple_mapped,$unmapped_mismatches,$unmapped_short,$unmapped_other" >> $output_file
done

#=============================================================
### step 4. Removing PCR duplication
#=============================================================
echo "Step 4: Deduplication"
inpath=${work_path}/s3_star_align

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
	{
		outpath=${work_path}/s4_bam_dedup/${sample}
		mkdir -p ${outpath}
		
		bam_in="${inpath}/${sample}/Aligned.sortedByCoord.out.bam"
		bai_in="${bam_in}.bai"
		bam_out="${outpath}/${sample}.dedup.bam"
		
		# ensure BAM index exists
		if [ ! -f "$bai_in" ] || [ "$bam_in" -nt "$bai_in" ]; then
			samtools index "$bam_in"
		fi
		
		# create a symbolic link since not using UMIs
		if [[ -e "$bam_out" ]]; then
			echo "Skipping symlink for $sample (already exists)"
		else
			ln -s "$bam_in" "$bam_out"
		fi
		# umi_tools dedup --paired -I ${inpath}/${sample}/Aligned.sortedByCoord.out.bam -S ${outpath}/${sample}.dedup.bam \
		# --multimapping-detection-method=NH --output-stats=${outpath}/deduplicated.txt --log=${outpath}/deduplication.log
	} &
	((job_count++))
	if ((job_count >= max_jobs)); then
		wait -n
		((job_count--))
	fi
done < $sample_file
wait

#=============================================================
### step 5. Bam file to Bed
#=============================================================
# Extract R2 reads that are successfully paired and aligned to the reference genome from the BAM file, then convert them to BED format.
echo "Step 5: Convert BAM to BED"
outdir=${work_path}/s5_R2_bed; mkdir -p $outdir
cd $outdir
bam_path=${work_path}/s4_bam_dedup

job_count=0
while IFS=$'\t' read -r sample fastq_1 fastq_2; do
    {
		final_bed=${outdir}/${sample}.dedup_R2.sorted.bed
		
		if [[ -f "$final_bed" ]]; then
			echo "Skipping $sample (BED already exists)"
			exit 0
		fi
		
		echo "Processing $sample..."
		samtools view -f 128 -F 4 -b ${bam_path}/${sample}/${sample}.dedup.bam > ${bam_path}/${sample}/${sample}.dedup_R2.bam
		samtools index ${bam_path}/${sample}/${sample}.dedup_R2.bam
		bedtools bamtobed -cigar -i ${bam_path}/${sample}/${sample}.dedup_R2.bam > ${sample}.dedup_R2.bed
		sort -k 1,1 ${sample}.dedup_R2.bed > ${sample}.dedup_R2.sorted.bed
		#gzip ${sample}.dedup_R2.sorted.bed
    }  &
	((job_count++))
	if ((job_count >= max_jobs)); then
		wait -n
		((job_count--))
	fi
done < $sample_file
wait

wc -l *.sorted.bed >> stats_deduped_reads.txt

# samtools view: -f 128, select read2 (second in pair); -F 4, Exclude mates that are unmapped
#=============================================================
### step 6. Defining LAPs (last aligned positions) and matching PASs
#=============================================================
# -misM 2: Allow up to 2 nt soft clipping at the ends of aligned reads.
# -dist 24: Assign LAPs to PASs in PolyA_DB within a ±24-nt window.
echo "Step 6: LAP/PAS analysis"
final_combined_file=${work_path}/s6_LAP/PAS_quant/cluster.all.reads.csv
if [[ -f "$final_combined_file" ]]; then
	echo "Skipping Step 6 (PAS/LAP output already exists)"
else
	indir=${work_path}/s5_R2_bed
	outdir=${work_path}/s6_LAP
	mkdir -p $outdir
	
	job_count=0
	while IFS=$'\t' read -r sample fastq_1 fastq_2; do
		{
			echo "Processing LAP/PAS for $sample..."
			Rscript ${script_path}/genome_LAPandPAS_define.R -bedLAP $indir/${sample}.dedup_R2.sorted.bed -out $outdir/${sample} -refPAS $refPAS_file -misM 2 -dist 24
		} &
		((job_count++))
		if ((job_count >= max_jobs)); then
			wait -n
			((job_count--))
		fi
	done < $sample_file
	wait 
	
	# combine the PAS quantification table
	mkdir -p ${work_path}/s6_LAP/PAS_quant
	cd ${work_path}/s6_LAP/PAS_quant
	mv ${work_path}/s6_LAP/*cluster.all.reads.csv ./
	Rscript ${script_path}/combine_all_sample_PAS_count_tables.R -csv ./ -out ./cluster.all.reads.csv --pattern "_cluster\\.all\\.reads\\.csv$"
	Rscript ${script_path}/combine_all_sample_PAS_count_tables.R -csv ${outdir} -out ./stats_PAS.csv --pattern ".*_stats\\.csv$"
fi
#=============================================================
### step 7. Generating bigwig files for PAS usage visualization
#=============================================================
echo "Step 7: Generate bigWig tracks"
indir=${work_path}/s6_LAP

for type in PASS_bw_LAP24 PASS_bw_LAP24_positon PASS_bw_LAP24_PAS mapped_read
do
outdir=${work_path}/s7_bigwig/${type}; mkdir -p $outdir; cd $outdir
echo "Working on ${type}..."

while IFS=$'\t' read -r sample fastq_1 fastq_2; do
  echo "Working on $sample..."
  ## Count total read number
  totalReadNum=`wc -l $indir/${sample}_${type}.bed | sed s/[[:blank:]].*//`
  echo "for file ${sample}_PASS_bw.bed, TotalReadNum=$totalReadNum"

  sort -k 1,1 $indir/${sample}_${type}.bed > $indir/${sample}_${type}.bed.sorted

  ## it is strand-specific 
  ## Generate bedgraph file
  echo "Generate bedgraph files for + and - strands..."
  genomeCoverageBed -bg -split -i $indir/${sample}_${type}.bed.sorted -strand '+' -g $chromsizes > $sample.plus.bedgraph
  genomeCoverageBed -bg -split -i $indir/${sample}_${type}.bed.sorted -strand '-' -g $chromsizes > $sample.minus.bedgraph
  
  #gzip $indir/${sample}_${type}.bed.sorted
  
  ## Normalize bedgraph counts
  echo "Normalize bedgraph counts..."
  ${script_path}/norm_bedgraph.pl -t $totalReadNum -i "$sample.plus.bedgraph"
  ${script_path}/norm_bedgraph.pl -t $totalReadNum -i "$sample.minus.bedgraph"

    ## give minus strand negative value
  awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2, $3, -$4}' $sample.minus.bedgraph.normalized > $sample.minus.bedgraph.normalized1 

  ## Convert to bigwig file
  echo "Convert to bigwig file..."
  bedGraphToBigWig $sample.plus.bedgraph.normalized  $chromsizes plus.$sample.$type.bw
  bedGraphToBigWig $sample.minus.bedgraph.normalized  $chromsizes minus.$sample.$type.bw
  bedGraphToBigWig $sample.minus.bedgraph.normalized1  $chromsizes mi.$sample.$type.bw

done < $sample_file
chmod 775 *.bw #change file permission
rm -f *.bedgraph
rm -f *.bedgraph.normalized
rm -f *.bedgraph.normalized1
done

echo "Pipeline complete!"
