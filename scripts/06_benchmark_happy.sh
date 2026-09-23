#!/usr/bin/env bash
#
# Step 06 — Benchmark each call set against the GIAB truth set with hap.py.
#
# hap.py (GA4GH benchmarking toolkit) performs haplotype-aware comparison via the
# vcfeval engine, so an indel represented differently but describing the same
# haplotype still counts as a true positive. It is run through Docker because the
# upstream image pins a working Python 2 / RTG toolchain that is painful to
# install natively.
#
# This step reproduces the original analysis: the comparison is confined to the
# GIAB high-confidence BED for the whole chromosome. That makes the denominator
# every confident variant on chr22, including the ~97% of the chromosome an exome
# capture never targets, so the recall reported here is a whole-chromosome recall
# and NOT on-target exome recall. See docs/LIMITATIONS.md, and run
# scripts/07_rerun_ontarget.sh for the target-restricted numbers.
#
# Results are written into results/happy/ and are small enough to commit.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd docker bcftools awk
require_file "${TRUTH_DIR}/${GIAB_TRUTH_VCF}" "run scripts/00_download_data.sh first"
require_file "${TRUTH_DIR}/${GIAB_TRUTH_BED}" "run scripts/00_download_data.sh first"

truth_vcf="${TRUTH_DIR}/${CHROM}_truth.vcf.gz"
conf_bed="${TRUTH_DIR}/${CHROM}_confident.bed"

# ------------------------------------------------------------------------------
# Restrict the truth set to the analysed chromosome
# ------------------------------------------------------------------------------

if [[ -s "$truth_vcf" ]]; then
  log "chr-restricted truth VCF exists — skipping"
else
  log "Extracting ${CHROM} from truth VCF"
  bcftools view -r "$CHROM" "${TRUTH_DIR}/${GIAB_TRUTH_VCF}" -Oz -o "$truth_vcf"
  bcftools index -t -f "$truth_vcf"
fi

if [[ -s "$conf_bed" ]]; then
  log "chr-restricted confident BED exists — skipping"
else
  log "Extracting ${CHROM} from confident-regions BED"
  awk -v c="$CHROM" '$1 == c' "${TRUTH_DIR}/${GIAB_TRUTH_BED}" > "$conf_bed"
  [[ -s "$conf_bed" ]] || die "no ${CHROM} intervals found in ${GIAB_TRUTH_BED}"
fi

# ------------------------------------------------------------------------------
# Run hap.py per depth
#
# The container sees two mounts: /data (working samples + truth) and /ref. Paths
# passed to hap.py must therefore be container paths, not host paths.
# ------------------------------------------------------------------------------

for depth in "${DEPTHS[@]}"; do
  base="$(bam_for_depth "$depth")"; base="${base%.bam}"
  name="$(basename "$base")"
  vcf="${base}.vcf.gz"

  require_file "$vcf" "run scripts/05_call_variants.sh first"

  log "Benchmarking ${depth}x with hap.py"
  docker run --rm \
    --platform "$HAPPY_PLATFORM" \
    -v "${SAMPLES_DIR}:/data" \
    -v "${TRUTH_DIR}:/truth" \
    -v "${REF_DIR}:/ref" \
    "$HAPPY_IMAGE" \
    /opt/hap.py/bin/hap.py \
      "/truth/$(basename "$truth_vcf")" \
      "/data/$(basename "$vcf")" \
      -f "/truth/$(basename "$conf_bed")" \
      -r "/ref/$(basename "$REFERENCE")" \
      -o "/data/${name}_happy" \
      --engine=vcfeval \
      --threads "$THREADS"

  cp "${base}_happy.summary.csv"  "${RESULTS_DIR}/happy/summary/${CHROM}_${depth}x.summary.csv"
  cp "${base}_happy.extended.csv" "${RESULTS_DIR}/happy/extended/${CHROM}_${depth}x.extended.csv"
  log "Done ${depth}x"
done

log "Step 06 complete. Tables in ${RESULTS_DIR}/happy/"
