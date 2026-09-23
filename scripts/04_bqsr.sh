#!/usr/bin/env bash
#
# Step 04 — Base Quality Score Recalibration.
#
# BQSR models the systematic error patterns a sequencer introduces into its
# reported base qualities. It does so by treating every mismatch that is NOT at a
# known variant site as a sequencing error, so the choice of known-sites resource
# directly determines what the model learns.
#
# IMPORTANT — this differs from the original submission on purpose. The original
# passed the GIAB benchmark VCF as --known-sites. That is the same VCF used to
# score the pipeline in step 06, which makes the recalibration circular: the true
# variant positions are excluded from the error model, and the model is tuned on
# knowledge the evaluation is supposed to test. This script uses dbSNP138 + Mills
# instead, per GATK Best Practices. Set BQSR_ALLOW_TRUTH_AS_KNOWN_SITES=1 to
# reproduce the original (circular) behaviour for comparison.
#
# Each depth gets a before/after covariate report and an AnalyzeCovariates PDF so
# the effect of recalibration is auditable.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd gatk samtools

known_args=()
if [[ "${BQSR_ALLOW_TRUTH_AS_KNOWN_SITES:-0}" == "1" ]]; then
  log "WARNING: using the benchmark truth set as BQSR known-sites (circular; original-submission behaviour)"
  known_args+=(--known-sites "${TRUTH_DIR}/${CHROM}_truth.vcf.gz")
else
  require_file "$KNOWN_SITES_DBSNP" "download the GATK hg38 resource bundle — see README"
  require_file "$KNOWN_SITES_MILLS" "download the GATK hg38 resource bundle — see README"
  known_args+=(--known-sites "$KNOWN_SITES_DBSNP" --known-sites "$KNOWN_SITES_MILLS")
fi

for depth in "${DEPTHS[@]}"; do
  base="$(bam_for_depth "$depth")"; base="${base%.bam}"
  in_bam="${base}.dedup.rg.bam"
  out_bam="${base}.bqsr.bam"
  before="${LOG_DIR}/$(basename "$base").recal_before.table"
  after="${LOG_DIR}/$(basename "$base").recal_after.table"
  plots="${LOG_DIR}/$(basename "$base").bqsr_covariates.pdf"

  require_file "$in_bam" "run scripts/03_preprocess.sh first"

  if [[ -s "$out_bam" ]]; then
    log "${depth}x already recalibrated — skipping"
    continue
  fi

  log "BQSR ${depth}x: building recalibration model"
  gatk --java-options "$GATK_JAVA_OPTS" BaseRecalibrator \
    -R "$REFERENCE" \
    -I "$in_bam" \
    -L "$CHROM" \
    "${known_args[@]}" \
    -O "$before"

  log "BQSR ${depth}x: applying recalibration"
  gatk --java-options "$GATK_JAVA_OPTS" ApplyBQSR \
    -R "$REFERENCE" \
    -I "$in_bam" \
    --bqsr-recal-file "$before" \
    --add-output-sam-program-record \
    --emit-original-quals \
    -O "$out_bam"

  log "BQSR ${depth}x: measuring residual bias"
  gatk --java-options "$GATK_JAVA_OPTS" BaseRecalibrator \
    -R "$REFERENCE" \
    -I "$out_bam" \
    -L "$CHROM" \
    "${known_args[@]}" \
    -O "$after"

  # AnalyzeCovariates needs R; a missing plotting stack should not abort a
  # six-hour run, so the failure is reported and the pipeline continues.
  gatk AnalyzeCovariates -before "$before" -after "$after" -plots "$plots" \
    || log "WARNING: AnalyzeCovariates failed for ${depth}x (R/ggplot2 missing?) — continuing"

  samtools index -@ "$THREADS" "$out_bam"
  log "Done ${depth}x -> ${out_bam}"
done

log "Step 04 complete."
