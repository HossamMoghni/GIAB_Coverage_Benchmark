#!/usr/bin/env bash
#
# Step 05 — Call variants with GATK HaplotypeCaller.
#
# HaplotypeCaller runs in GVCF mode and GenotypeGVCFs then assigns genotypes.
# The two-stage route is used rather than direct VCF output because the GVCF
# records reference confidence at non-variant sites, which keeps the per-depth
# calls comparable: a site absent from a shallow call set is distinguishable from
# a site that was never examined.
#
# --pcr-indel-model NONE is correct here because the library is PCR-free.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd gatk bgzip bcftools

for depth in "${DEPTHS[@]}"; do
  base="$(bam_for_depth "$depth")"; base="${base%.bam}"
  in_bam="${base}.bqsr.bam"
  gvcf="${base}.g.vcf.gz"
  vcf="${base}.vcf.gz"

  require_file "$in_bam" "run scripts/04_bqsr.sh first"

  if [[ -s "$vcf" ]]; then
    log "${depth}x already called — skipping"
    continue
  fi

  log "HaplotypeCaller ${depth}x (GVCF mode)"
  gatk --java-options "$GATK_JAVA_OPTS" HaplotypeCaller \
    -R "$REFERENCE" \
    -I "$in_bam" \
    -L "$CHROM" \
    --pcr-indel-model NONE \
    --emit-ref-confidence GVCF \
    -O "$gvcf"

  log "GenotypeGVCFs ${depth}x"
  gatk --java-options "$GATK_JAVA_OPTS" GenotypeGVCFs \
    -R "$REFERENCE" \
    -V "$gvcf" \
    -L "$CHROM" \
    --max-alternate-alleles 6 \
    -O "$vcf"

  bcftools index -t -f "$vcf"
  log "Done ${depth}x -> ${vcf}"
done

log "Step 05 complete."
