#!/bin/bash
# Variant Calling Benchmarking Pipeline
# Project: Evaluating Variant Calling Performance at Different Depths of Coverage
# Student Name: Hossam Hatem

# ==============================================================================
# SECTION 1: DATA DOWNLOAD
# ==============================================================================

echo "Step 1: Downloading HG002 Exome Data"

# Download BAM file from GIAB
# Source: https://github.com/genome-in-a-bottle/giab_data_indexes
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam

# Download BAM index
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai

echo ""
echo "Step 2: Downloading GRCh38 Reference"

# Download GRCh38 reference genome
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Decompress and rename
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa

# ==============================================================================
# SECTION 2: BAM TO FASTQ CONVERSION & REALIGNMENT
# ==============================================================================

echo ""
echo "Step 1: Sort BAM by Read Name"

samtools sort -n -@ 8 \
  -o full_HG002.sorted.bam \
  151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam

echo ""
echo "Step 2: Convert BAM to FASTQ"

samtools fastq -@ 8 \
  -1 HG002_exome_full_R1.fastq.gz \
  -2 HG002_exome_full_R2.fastq.gz \
  -0 /dev/null \
  -s /dev/null \
  -n full_HG002.sorted.bam

echo ""
echo "Step 3: Index Reference for BWA"
echo ""

bwa index GRCh38.fa

echo ""
echo "Step 4: Align FASTQ to GRCh38"
echo ""

bwa mem -t 8 \
  -R "@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
  "/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/ref/GRCh38.fa" \
  HG002_exome_full_R1.fastq.gz \
  HG002_exome_full_R2.fastq.gz | \
samtools view -@ 8 -b -o full_38_aligned.bam -

echo ""
echo "Step 5: Sort and Index BAM"
echo ""

samtools sort -@ 8 -o aligned_sorted.bam full_38_aligned.bam
samtools index aligned_sorted.bam

# Verify chromosome names
echo ""
echo ">>> Checking chromosome nomenclature:"
samtools view -H aligned_sorted.bam | grep '@SQ' | head -5

# ==============================================================================
# SECTION 3: EXTRACT CHROMOSOME 22 & COVERAGE SUBSAMPLING
# ==============================================================================

echo ""
echo ""
echo "Step 1: Extract Chromosome 22"
echo ""

samtools view -b -@ 8 \
  aligned_sorted.bam \
  chr22 > HG002_exome_chr22.bam

echo ""
echo "Step 2: Index chr22 BAM"
echo ""

samtools index HG002_exome_chr22.bam

echo ""
echo ""
echo "Step 3: Calculate Coverage Statistics"

# Overall average coverage
echo ">>> Overall average coverage:"
samtools depth HG002_exome_chr22.bam | \
  awk '{sum+=$3} END { print "Average coverage:", sum/NR }'

# Coverage for regions ≥5x
echo ""
echo ">>> Coverage for regions ≥5x:"
samtools depth HG002_exome_chr22.bam | \
  awk '$3 >= 5 {sum+=$3; count++} END {
    print "Average coverage (≥5x):", sum/count
    print "Positions with ≥5x:", count
  }'

echo ""
echo "Step 4: Subsample to Different Coverages"

# Subsample to 40x coverage (65.3% of reads)
echo ">>> Subsampling to 40x coverage..."
samtools view -b -s 42.653 -@ 8 \
  HG002_exome_chr22.bam > HG002_exome_chr22_40x.bam

# Subsample to 10x coverage (16.3% of reads)
echo ">>> Subsampling to 10x coverage..."
samtools view -b -s 42.163 -@ 8 \
  HG002_exome_chr22.bam > HG002_exome_chr22_10x.bam

# Subsample to 2x coverage (3.3% of reads)
echo ">>> Subsampling to 2x coverage..."
samtools view -b -s 42.033 -@ 8 \
  HG002_exome_chr22.bam > HG002_exome_chr22_2x.bam

# Index all subsampled BAMs
echo ">>> Indexing subsampled BAMs..."
samtools index HG002_exome_chr22_40x.bam
samtools index HG002_exome_chr22_10x.bam
samtools index HG002_exome_chr22_2x.bam

echo "Subsampling complete!"

# ==============================================================================
# SECTION 4: MARK DUPLICATES & ADD READ GROUPS
# ==============================================================================

echo ""
echo ""
echo "SECTION 4: MARK DUPLICATES & ADD READ GROUPS"
echo ""

echo ""
echo ">>> Step 4.1: Mark Duplicates"

# Loop through all coverage samples
for sample in HG002_exome_chr22_*x.bam; do
  name=${sample%.bam}
  depth=$(echo $name | grep -o '[0-9]\+x')
  
  echo ">>> Marking duplicates in ${depth}..."
  
  picard MarkDuplicates \
    INPUT=$sample \
    OUTPUT=${name}.dedup.bam \
    METRICS_FILE=${name}.metrics.txt \
    VALIDATION_STRINGENCY=SILENT
  
  # Index the deduplicated BAM
  samtools index ${name}.dedup.bam
  
  echo "Done with ${depth}"
done

echo ""
echo ">>> Step 4.2: Add Read Groups"

# Add read groups to all deduplicated BAMs
for sample in HG002_exome_chr22_*.dedup.bam; do
  name=${sample%.dedup.bam}
  depth=$(echo $name | grep -o '[0-9]\+x')
  
  echo ">>> Adding read groups to ${depth}..."
  
  picard AddOrReplaceReadGroups \
    INPUT=$sample \
    OUTPUT=${name}.dedup.rg.bam \
    RGID=${depth} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=HG002 \
    VALIDATION_STRINGENCY=SILENT
  
  samtools index ${name}.dedup.rg.bam
  echo "Done with ${depth}"
done


# ==============================================================================
# SECTION 5: PREPARE TRUTH SET & BQSR
# ==============================================================================

echo ""
echo "SECTION 5: PREPARE TRUTH SET & BQSR"

echo ""
echo ">>> Step 5.1: Prepare GIAB Truth Set"

# Remove old index if exists
rm -f HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

# Create new index
bcftools index -t HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

echo ""
echo ">>> Step 5.2: Extract chr22 from Truth Set"

# Extract chr22
bcftools view -r chr22 \
  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  -Oz -o HG002_chr22_truth.vcf.gz

# Index the chr22 file
bcftools index -t HG002_chr22_truth.vcf.gz

echo ""
echo ">>> Step 5.3: Extract chr22 Confident Regions"

# Extract chr22 from confident regions BED
awk '$1=="chr22" {print}' \
  HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  > HG002_chr22_confident.bed

echo ""
echo ">>> Step 5.4: Base Quality Score Recalibration (BQSR)"

# BQSR loop for all samples
for sample in HG002_exome_chr22_*.dedup.rg.bam; do
  name=${sample%.dedup.rg.bam}
  depth=$(echo $name | grep -o '[0-9]\+x')
  
  echo ""
  echo ">>> Processing ${depth} for BQSR..."
  
  # First pass: BaseRecalibrator (before BQSR)
  echo "   - Running BaseRecalibrator (before)..."
  gatk --java-options "-Xmx10G" BaseRecalibrator \
    -R "/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/ref/GRCh38.fa" \
    -I $sample \
    --known-sites HG002_chr22_truth.vcf.gz \
    -O ${name}.report
  
  # Apply BQSR
  echo "   - Applying BQSR..."
  gatk --java-options "-Xmx10G" ApplyBQSR \
    -R "/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/ref/GRCh38.fa" \
    -I $sample \
    --bqsr-recal-file ${name}.report \
    -O ${name}.bqsr.bam \
    --add-output-sam-program-record \
    --emit-original-quals
  
  # Second pass: BaseRecalibrator (after BQSR)
  echo "   - Running BaseRecalibrator (after)..."
  gatk --java-options "-Xmx10G" BaseRecalibrator \
    -R "/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/ref/GRCh38.fa" \
    -I ${name}.bqsr.bam \
    --known-sites HG002_chr22_truth.vcf.gz \
    -O ${name}.report2
  
  # Generate BQSR comparison plots
  gatk AnalyzeCovariates \
    -before ${name}.report \
    -after ${name}.report2 \
    -plots ${name}.bqsr_plots.pdf
  
  # Index the BQSR BAM
  samtools index ${name}.bqsr.bam
  
  echo "Done with $name"
done


Step 1: Variant Calling (HaplotypeCaller + GenotypeGVCFs)
for sample in *.bqsr.bam; do
  name=${sample%.bqsr.bam}
  
  echo "Calling variants for $name"
  
  # HaplotypeCaller (GVCF mode) - chr22 only
  gatk --java-options "-Xmx10G" HaplotypeCaller \
    -R "/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/ref/GRCh38.fa" \
    -I "$sample" \
    -L chr22 \
    --pcr-indel-model NONE \
    --emit-ref-confidence GVCF \
    -O ${name}.g.vcf
  
  # GenotypeGVCFs (convert GVCF to VCF) - chr22 only
  echo "   - Running GenotypeGVCFs..."
  gatk --java-options "-Xmx10G" GenotypeGVCFs \
    -R "/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/ref/GRCh38.fa" \
    -V ${name}.g.vcf \
    -L chr22 \
    --max-alternate-alleles 6 \
    -O ${name}.vcf
  
  echo "Done with ${depth}"
done

echo ""
echo ">>> Compressing and Indexing VCF Files"

# Compress and index VCF files for hap.py
for sample in HG002_exome_chr22_*.vcf; do
  if [ -f "$sample" ]; then
    echo ">>> Compressing $sample..."
    bgzip -f $sample
    bcftools index -t -f ${sample}.gz
  fi
done


# ==============================================================================
# SECTION 7: BENCHMARKING WITH HAP.PY
# ==============================================================================

echo ""
echo "SECTION 7: BENCHMARKING WITH HAP.PY"

# Set directories
SAMPLES_DIR="/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/samples"
REF_DIR="/Users/hossam/Documents/NU Master's/3rd Semester (Fall 2025)/Advanced NGS/Project/ref"

echo ""
echo ">>> Running hap.py benchmarking via Docker"

# Loop through all samples
for sample in HG002_exome_chr22_*.vcf.gz; do
  name=${sample%.vcf.gz}
  depth=$(echo $name | grep -o '[0-9]\+x')
  
  echo ""
  echo ">>> Benchmarking ${depth}..."
  
  docker run --rm \
    --platform linux/amd64 \
    -v "${SAMPLES_DIR}:/data" \
    -v "${REF_DIR}:/ref" \
    jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
    /data/HG002_chr22_truth.vcf.gz \
    /data/${sample} \
    -f /data/HG002_chr22_confident.bed \
    -r /ref/GRCh38.fa \
    -o /data/${name}_happy \
    --engine=vcfeval \
    --threads 4
  
  echo "Done with ${depth}"
done