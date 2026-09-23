#!/usr/bin/env bash
#
# Step 02 — Extract chromosome 22 and build the coverage series.
#
# The full chr22 extraction is the deepest point in the series (measured mean
# depth over covered positions is ~61x, reported as 60x). Lower depths are drawn
# from it with `samtools view -s`, using a fixed seed so the read subsets are
# reproducible across runs.
#
# Unlike the original submission, the subsampling fractions are NOT hard-coded:
# they are computed from the measured full depth, so the series stays correct if
# the input data or the target depths change. The full-depth BAM is also copied
# into the series under its own depth label, so downstream steps process all four
# depths uniformly.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd samtools awk
sorted="${PROJECT_DIR}/HG002_exome_GRCh38.sorted.bam"
require_file "$sorted" "run scripts/01_realign_grch38.sh first"

full_bam="$(bam_for_depth "$FULL_DEPTH")"
depth_stats="${LOG_DIR}/chr22_depth_stats.txt"

# ------------------------------------------------------------------------------
# Extract chr22
# ------------------------------------------------------------------------------

if [[ -s "$full_bam" ]]; then
  log "chr22 BAM exists — skipping extraction"
else
  log "Extracting ${CHROM}"
  samtools view -b -@ "$THREADS" "$sorted" "$CHROM" > "$full_bam"
  samtools index -@ "$THREADS" "$full_bam"
fi

# ------------------------------------------------------------------------------
# Measure depth
#
# Two numbers are reported, because they answer different questions:
#   * mean over covered positions only  — the depth an exome read "feels"
#   * mean over positions at >= 5x      — excludes off-target single-read noise
#
# The first is used to derive subsampling fractions, matching how the original
# analysis defined its depth targets.
# ------------------------------------------------------------------------------

log "Measuring coverage on ${CHROM} (this reads the whole BAM)"

samtools depth "$full_bam" | awk '
  { sum += $3; n++ }
  $3 >= 5 { sum5 += $3; n5++ }
  END {
    printf "mean_depth_covered\t%.4f\n", (n  ? sum  / n  : 0)
    printf "covered_positions\t%d\n",    n
    printf "mean_depth_ge5x\t%.4f\n",    (n5 ? sum5 / n5 : 0)
    printf "positions_ge5x\t%d\n",       n5
  }' > "$depth_stats"

cat "$depth_stats" >&2

full_depth_measured=$(awk '$1=="mean_depth_covered" {print $2}' "$depth_stats")
[[ -n "$full_depth_measured" ]] || die "failed to measure depth from ${full_bam}"

log "Measured full-coverage depth: ${full_depth_measured}x (labelled ${FULL_DEPTH}x)"

# ------------------------------------------------------------------------------
# Subsample to each lower depth
# ------------------------------------------------------------------------------

for depth in "${DEPTHS[@]}"; do
  [[ "$depth" == "$FULL_DEPTH" ]] && continue

  out="$(bam_for_depth "$depth")"
  if [[ -s "$out" ]]; then
    log "${depth}x BAM exists — skipping"
    continue
  fi

  # samtools -s takes SEED.FRACTION as a single float, e.g. 42.653 means
  # "seed 42, keep 65.3% of reads". Build that string from the measured depth.
  fraction=$(awk -v want="$depth" -v have="$full_depth_measured" \
    'BEGIN { f = want / have; if (f > 0.999) f = 0.999; printf "%.3f", f }')

  # A requested depth at or above the measured full depth cannot be produced by
  # downsampling. Fail loudly rather than silently emitting a near-full BAM
  # mislabelled as something shallower.
  awk -v want="$depth" -v have="$full_depth_measured" \
    'BEGIN { exit (want < have) ? 0 : 1 }' \
    || die "requested ${depth}x exceeds measured full depth ${full_depth_measured}x — adjust DEPTHS or FULL_DEPTH"

  log "Subsampling to ${depth}x (keeping ${fraction} of reads)"
  samtools view -b -@ "$THREADS" \
    -s "${SUBSAMPLE_SEED}${fraction#0}" \
    "$full_bam" > "$out"
  samtools index -@ "$THREADS" "$out"
done

log "Step 02 complete. Coverage series in ${SAMPLES_DIR}:"
ls -lh "${SAMPLES_DIR}"/*.bam >&2
