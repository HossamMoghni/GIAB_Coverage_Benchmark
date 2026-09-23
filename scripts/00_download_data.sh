#!/usr/bin/env bash
#
# Step 00 — Download input data.
#
#   * GIAB HG002 whole-exome BAM (GRCh37-aligned) + index
#   * GRCh38 no-alt analysis-set reference genome
#   * GIAB v4.2.1 benchmark truth VCF + high-confidence BED
#
# All downloads are resumable (`wget -c` / `curl -C -`), so re-running this
# script after an interrupted transfer picks up where it stopped. Files that are
# already complete are skipped.
#
# Disk required: ~60 GB.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

require_cmd wget curl gunzip

# ------------------------------------------------------------------------------
# HG002 exome BAM
# ------------------------------------------------------------------------------

log "Downloading HG002 exome BAM and index into ${RAW_DIR}"

wget -c -P "$RAW_DIR" "${GIAB_EXOME_BASE}/${GIAB_EXOME_BAM}"
wget -c -P "$RAW_DIR" "${GIAB_EXOME_BASE}/${GIAB_EXOME_BAM%.bam}.bai"

# ------------------------------------------------------------------------------
# GRCh38 reference
# ------------------------------------------------------------------------------

if [[ -s "$REFERENCE" ]]; then
  log "Reference already present: ${REFERENCE} — skipping"
else
  log "Downloading GRCh38 reference"
  ref_gz="${REF_DIR}/$(basename "$REFERENCE_URL")"
  curl -L -C - -o "$ref_gz" "$REFERENCE_URL"

  log "Decompressing reference"
  gunzip -c "$ref_gz" > "$REFERENCE"
  rm -f "$ref_gz"
fi

# ------------------------------------------------------------------------------
# GIAB truth set
# ------------------------------------------------------------------------------

log "Downloading GIAB v4.2.1 truth set"

curl -L -C - -o "${TRUTH_DIR}/${GIAB_TRUTH_VCF}"     "${GIAB_TRUTH_BASE}/${GIAB_TRUTH_VCF}"
curl -L -C - -o "${TRUTH_DIR}/${GIAB_TRUTH_VCF}.tbi" "${GIAB_TRUTH_BASE}/${GIAB_TRUTH_VCF}.tbi"
curl -L -C - -o "${TRUTH_DIR}/${GIAB_TRUTH_BED}"     "${GIAB_TRUTH_BASE}/${GIAB_TRUTH_BED}"

log "Step 00 complete. Raw data in ${RAW_DIR}, reference in ${REF_DIR}, truth set in ${TRUTH_DIR}"
