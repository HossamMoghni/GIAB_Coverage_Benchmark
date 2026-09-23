#!/usr/bin/env bash
#
# Step 03 — Mark duplicates and normalise read groups.
#
# Duplicates are marked per depth, AFTER subsampling. That ordering matters: a
# downsampled library has genuinely fewer duplicate pairs, so marking before the
# split would carry the full-depth duplicate rate into every shallow sample and
# understate their usable coverage.
#
# Read groups are then rewritten so each depth carries a distinct RGID while
# sharing one sample name — GATK requires read groups, and the per-depth ID makes
# the provenance of any downstream BAM obvious.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd samtools picard

for depth in "${DEPTHS[@]}"; do
  in_bam="$(bam_for_depth "$depth")"
  base="${in_bam%.bam}"
  dedup="${base}.dedup.bam"
  final="${base}.dedup.rg.bam"

  require_file "$in_bam" "run scripts/02_subsample_coverage.sh first"

  if [[ -s "$final" ]]; then
    log "${depth}x already preprocessed — skipping"
    continue
  fi

  log "Marking duplicates at ${depth}x"
  picard MarkDuplicates \
    INPUT="$in_bam" \
    OUTPUT="$dedup" \
    METRICS_FILE="${LOG_DIR}/$(basename "$base").duplicate_metrics.txt" \
    VALIDATION_STRINGENCY=SILENT

  log "Adding read groups at ${depth}x"
  picard AddOrReplaceReadGroups \
    INPUT="$dedup" \
    OUTPUT="$final" \
    RGID="${depth}x" \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM="$SAMPLE_ID" \
    VALIDATION_STRINGENCY=SILENT

  samtools index -@ "$THREADS" "$final"
  rm -f "$dedup"

  log "Done ${depth}x -> ${final}"
done

log "Step 03 complete."
