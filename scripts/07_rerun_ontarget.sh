#!/usr/bin/env bash
#
# Step 07 — Re-benchmark restricted to the exome capture target (correction).
#
# Why this exists
# ---------------
# Step 06 scores the call sets against every GIAB high-confidence variant on
# chr22. Chromosome 22 has roughly 39 Mb of high-confidence sequence, but an
# Agilent SureSelect v5 capture targets only about 1 Mb of it. Variants outside
# the capture are unreachable by this library at any depth, yet step 06 counts
# every one of them as a false negative. The result is a recall ceiling of very
# roughly 3%, which is why the headline numbers look catastrophic (25% SNP recall
# at 60x) and why the reported precision is far more informative than the recall.
#
# Passing the capture BED to hap.py as --target-regions (-T) makes the
# denominator "confident variants this assay could actually see", which is the
# quantity the coverage-versus-accuracy question is really about. -T restricts
# scoring; -f (the confident BED) still defines which calls are counted at all.
#
# Requirements
# ------------
# An exome target BED in GRCh38 coordinates, restricted to the analysed
# chromosome, at $EXOME_TARGET_BED. For this library that is the Agilent
# SureSelect Human All Exon v5 "Covered" BED, obtained from Agilent SureDesign
# (free account required); lift it to GRCh38 if the download is GRCh37.
#
# Outputs land in results/happy_ontarget/ and are written alongside — not over —
# the whole-chromosome tables, so the two views stay comparable.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd docker awk
require_file "$EXOME_TARGET_BED" "exome capture BED in GRCh38 coordinates — see the header of this script"

truth_vcf="${TRUTH_DIR}/${CHROM}_truth.vcf.gz"
conf_bed="${TRUTH_DIR}/${CHROM}_confident.bed"
require_file "$truth_vcf" "run scripts/06_benchmark_happy.sh first"
require_file "$conf_bed"  "run scripts/06_benchmark_happy.sh first"

# Keep only the analysed chromosome, in case a genome-wide BED was supplied.
target_bed="${TRUTH_DIR}/${CHROM}_exome_target.bed"
awk -v c="$CHROM" '$1 == c' "$EXOME_TARGET_BED" > "$target_bed"
[[ -s "$target_bed" ]] || die "no ${CHROM} intervals in ${EXOME_TARGET_BED} — is it GRCh38 with UCSC-style contig names?"

out_summary="${RESULTS_DIR}/happy_ontarget/summary"
out_extended="${RESULTS_DIR}/happy_ontarget/extended"
mkdir -p "$out_summary" "$out_extended"

log "On-target BED: $(wc -l < "$target_bed") intervals covering $(awk '{n += $3 - $2} END {print n}' "$target_bed") bp"

for depth in "${DEPTHS[@]}"; do
  base="$(bam_for_depth "$depth")"; base="${base%.bam}"
  name="$(basename "$base")"
  vcf="${base}.vcf.gz"

  require_file "$vcf" "run scripts/05_call_variants.sh first"

  log "On-target benchmarking ${depth}x"
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
      -T "/truth/$(basename "$target_bed")" \
      -r "/ref/$(basename "$REFERENCE")" \
      -o "/data/${name}_happy_ontarget" \
      --engine=vcfeval \
      --threads "$THREADS"

  cp "${base}_happy_ontarget.summary.csv"  "${out_summary}/${CHROM}_${depth}x.summary.csv"
  cp "${base}_happy_ontarget.extended.csv" "${out_extended}/${CHROM}_${depth}x.extended.csv"
  log "Done ${depth}x"
done

log "Step 07 complete. On-target tables in ${RESULTS_DIR}/happy_ontarget/"
log "Regenerate figures with: python analysis/plot_benchmarks.py --results-dir ${RESULTS_DIR}/happy_ontarget"
