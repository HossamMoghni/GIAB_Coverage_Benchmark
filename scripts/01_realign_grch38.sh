#!/usr/bin/env bash
#
# Step 01 — Realign the GRCh37 exome BAM to GRCh38.
#
# The GIAB exome BAM is aligned to GRCh37, but the v4.2.1 benchmark truth set
# used in step 06 is in GRCh38 coordinates. Rather than lifting coordinates over
# (which loses reads in rearranged regions), the reads are extracted back to
# FASTQ and realigned from scratch with BWA-MEM.
#
#   name-sorted BAM -> paired FASTQ -> bwa mem GRCh38 -> coordinate-sorted BAM
#
# Runtime: several hours on 8 threads. Disk required: ~150 GB peak.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd samtools bwa
require_file "${RAW_DIR}/${GIAB_EXOME_BAM}" "run scripts/00_download_data.sh first"
require_file "$REFERENCE" "run scripts/00_download_data.sh first"

namesorted="${PROJECT_DIR}/HG002_exome.namesorted.bam"
fq1="${PROJECT_DIR}/HG002_exome_R1.fastq.gz"
fq2="${PROJECT_DIR}/HG002_exome_R2.fastq.gz"
aligned="${PROJECT_DIR}/HG002_exome_GRCh38.bam"
sorted="${PROJECT_DIR}/HG002_exome_GRCh38.sorted.bam"

# ------------------------------------------------------------------------------
# Name-sort, so that mates stay adjacent for `samtools fastq`
# ------------------------------------------------------------------------------

if [[ -s "$namesorted" ]]; then
  log "Name-sorted BAM exists — skipping"
else
  log "Name-sorting source BAM"
  samtools sort -n -@ "$THREADS" \
    -o "$namesorted" \
    "${RAW_DIR}/${GIAB_EXOME_BAM}"
fi

# ------------------------------------------------------------------------------
# BAM -> paired FASTQ
# ------------------------------------------------------------------------------

if [[ -s "$fq1" && -s "$fq2" ]]; then
  log "FASTQ pair exists — skipping conversion"
else
  log "Converting BAM to paired FASTQ"
  samtools fastq -@ "$THREADS" \
    -1 "$fq1" \
    -2 "$fq2" \
    -0 /dev/null \
    -s /dev/null \
    -n "$namesorted"
fi

# ------------------------------------------------------------------------------
# Index the reference for BWA (one-off, ~1 hour)
# ------------------------------------------------------------------------------

if [[ -s "${REFERENCE}.bwt" ]]; then
  log "BWA index exists — skipping"
else
  log "Building BWA index (this takes roughly an hour)"
  bwa index "$REFERENCE"
fi

# GATK also needs a .fai and a sequence dictionary downstream.
[[ -s "${REFERENCE}.fai" ]]          || samtools faidx "$REFERENCE"
[[ -s "${REFERENCE%.fa}.dict" ]]     || samtools dict "$REFERENCE" -o "${REFERENCE%.fa}.dict"

# ------------------------------------------------------------------------------
# Align to GRCh38
# ------------------------------------------------------------------------------

if [[ -s "$sorted" ]]; then
  log "Coordinate-sorted GRCh38 BAM exists — skipping alignment"
else
  log "Aligning to GRCh38 with BWA-MEM"
  bwa mem -t "$THREADS" \
    -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "$REFERENCE" \
    "$fq1" "$fq2" \
  | samtools view -@ "$THREADS" -b -o "$aligned" -

  log "Coordinate-sorting and indexing"
  samtools sort -@ "$THREADS" -o "$sorted" "$aligned"
  samtools index -@ "$THREADS" "$sorted"
  rm -f "$aligned"
fi

log "Chromosome nomenclature check (expect UCSC-style 'chrN'):"
samtools view -H "$sorted" | grep '@SQ' | head -5 >&2

log "Step 01 complete: ${sorted}"
