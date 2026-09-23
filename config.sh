#!/usr/bin/env bash
#
# Central configuration for the variant-calling benchmarking pipeline.
#
# Every script in scripts/ sources this file. Override any value by exporting it
# before you run the pipeline, e.g.:
#
#     THREADS=16 PROJECT_DIR=/scratch/ngs ./scripts/run_all.sh
#
# Nothing here writes to disk; it only defines variables and creates the
# directory layout the rest of the pipeline expects.

set -euo pipefail

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

# Root of this repository (resolved from the location of this file, so the
# pipeline works no matter where it is invoked from).
REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Working directory for large intermediate data. This is deliberately OUTSIDE
# the repository by default: BAM, FASTQ and reference files are tens of
# gigabytes and must never be committed to git.
PROJECT_DIR="${PROJECT_DIR:-${HOME}/ngs_benchmarking_work}"

REF_DIR="${REF_DIR:-${PROJECT_DIR}/ref}"
RAW_DIR="${RAW_DIR:-${PROJECT_DIR}/raw}"
SAMPLES_DIR="${SAMPLES_DIR:-${PROJECT_DIR}/samples}"
TRUTH_DIR="${TRUTH_DIR:-${PROJECT_DIR}/truth}"
LOG_DIR="${LOG_DIR:-${PROJECT_DIR}/logs}"

# Final benchmarking tables are small and DO belong in the repository.
RESULTS_DIR="${RESULTS_DIR:-${REPO_DIR}/results}"

# ------------------------------------------------------------------------------
# Reference genome
# ------------------------------------------------------------------------------

REFERENCE="${REFERENCE:-${REF_DIR}/GRCh38.fa}"
REFERENCE_URL="${REFERENCE_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz}"

# ------------------------------------------------------------------------------
# Input data: GIAB HG002 whole-exome BAM (aligned to GRCh37, realigned in step 01)
# ------------------------------------------------------------------------------

GIAB_EXOME_BASE="${GIAB_EXOME_BASE:-ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome}"
GIAB_EXOME_BAM="${GIAB_EXOME_BAM:-151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam}"

# GIAB v4.2.1 benchmark truth set (VCF + high-confidence BED).
GIAB_TRUTH_BASE="${GIAB_TRUTH_BASE:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38}"
GIAB_TRUTH_VCF="${GIAB_TRUTH_VCF:-HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
GIAB_TRUTH_BED="${GIAB_TRUTH_BED:-HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"

# Known-sites resources for BQSR (step 04).
#
# These must be INDEPENDENT of the benchmark truth set. The original submission
# passed the GIAB truth VCF here, which masks exactly the positions the pipeline
# is later scored on and makes the recalibration circular. dbSNP + Mills are the
# GATK Best Practices resources for this step.
KNOWN_SITES_DBSNP="${KNOWN_SITES_DBSNP:-${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf}"
KNOWN_SITES_MILLS="${KNOWN_SITES_MILLS:-${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}"
GATK_BUNDLE_URL="${GATK_BUNDLE_URL:-https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0}"

# Exome capture target intervals (Agilent SureSelect v5 for this library).
# Used by scripts/07_rerun_ontarget.sh to restrict benchmarking to captured
# regions. See docs/LIMITATIONS.md for why this matters.
EXOME_TARGET_BED="${EXOME_TARGET_BED:-${REF_DIR}/agilent_sureselect_v5_chr22.bed}"

# ------------------------------------------------------------------------------
# Analysis parameters
# ------------------------------------------------------------------------------

CHROM="${CHROM:-chr22}"
SAMPLE_ID="${SAMPLE_ID:-HG002}"
THREADS="${THREADS:-8}"
GATK_JAVA_OPTS="${GATK_JAVA_OPTS:--Xmx10G}"

# Fixed seed for samtools subsampling, so every run reproduces the same reads.
SUBSAMPLE_SEED="${SUBSAMPLE_SEED:-42}"

# Target depths in x. "60" is the full chr22 extraction (measured mean depth on
# covered positions is ~61x), so it is not subsampled; the others are drawn from
# it. Fractions are recomputed at run time by scripts/02_subsample_coverage.sh
# from the measured full depth rather than being hard-coded.
if [[ -z "${DEPTHS+set}" ]]; then
  DEPTHS=(2 10 40 60)
fi

# The depth label that corresponds to the un-subsampled full chr22 extraction.
FULL_DEPTH="${FULL_DEPTH:-60}"

# hap.py Docker image (GA4GH benchmarking toolkit).
HAPPY_IMAGE="${HAPPY_IMAGE:-jmcdani20/hap.py:v0.3.12}"
HAPPY_PLATFORM="${HAPPY_PLATFORM:-linux/amd64}"

# ------------------------------------------------------------------------------
# Helpers shared by all steps
# ------------------------------------------------------------------------------

log()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()  { printf '[%s] ERROR: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; exit 1; }

# require_cmd tool [tool...] — fail early with a readable message instead of
# dying halfway through a 6-hour run.
require_cmd() {
  local missing=()
  for cmd in "$@"; do
    command -v "$cmd" >/dev/null 2>&1 || missing+=("$cmd")
  done
  if [[ ${#missing[@]} -gt 0 ]]; then
    die "missing required tool(s): ${missing[*]} — see environment.yml"
  fi
}

# require_file path [description] — verify an input exists before using it.
require_file() {
  [[ -s "$1" ]] || die "missing or empty input: $1${2:+ ($2)}"
}

# bam_for_depth 40 -> /path/samples/HG002_chr22_40x.bam
bam_for_depth() { printf '%s/%s_%s_%sx.bam' "$SAMPLES_DIR" "$SAMPLE_ID" "$CHROM" "$1"; }

mkdir -p "$REF_DIR" "$RAW_DIR" "$SAMPLES_DIR" "$TRUTH_DIR" "$LOG_DIR" \
         "$RESULTS_DIR/happy/summary" "$RESULTS_DIR/happy/extended" "$RESULTS_DIR/figures"
