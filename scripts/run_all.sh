#!/usr/bin/env bash
#
# Run the whole pipeline end to end.
#
#   ./scripts/run_all.sh              # steps 00-06, then figures
#   ./scripts/run_all.sh 02 03 04     # only the named steps
#
# Every step is resumable: completed outputs are detected and skipped, so an
# interrupted run can simply be restarted. Step 07 (on-target correction) is not
# part of the default run because it needs an exome capture BED you must obtain
# separately — see the header of scripts/07_rerun_ontarget.sh.
#
# Expect roughly 12-24 hours on 8 cores and ~200 GB of scratch disk.

# shellcheck source=../config.sh
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ALL_STEPS=(
  "00_download_data.sh"
  "01_realign_grch38.sh"
  "02_subsample_coverage.sh"
  "03_preprocess.sh"
  "04_bqsr.sh"
  "05_call_variants.sh"
  "06_benchmark_happy.sh"
)

# Select steps: either the ones named on the command line, or all of them.
steps=()
if [[ $# -gt 0 ]]; then
  for want in "$@"; do
    matched=""
    for step in "${ALL_STEPS[@]}" "07_rerun_ontarget.sh"; do
      if [[ "$step" == "$want" || "$step" == "${want}_"* ]]; then
        steps+=("$step"); matched=1; break
      fi
    done
    [[ -n "$matched" ]] || die "unknown step: ${want} (expected one of: ${ALL_STEPS[*]})"
  done
else
  steps=("${ALL_STEPS[@]}")
fi

started=$(date +%s)

for step in "${steps[@]}"; do
  log "================ ${step} ================"
  step_started=$(date +%s)

  bash "${SCRIPT_DIR}/${step}" 2>&1 | tee -a "${LOG_DIR}/${step%.sh}.log"

  # tee sits at the end of the pipe, so check the step's own status.
  status=${PIPESTATUS[0]}
  [[ $status -eq 0 ]] || die "${step} failed with exit status ${status} — see ${LOG_DIR}/${step%.sh}.log"

  log "${step} finished in $(( ($(date +%s) - step_started) / 60 )) min"
done

# Figures are cheap to rebuild and depend only on the committed CSVs.
if command -v python3 >/dev/null 2>&1; then
  log "================ figures ================"
  python3 "${REPO_DIR}/analysis/plot_benchmarks.py" || log "WARNING: figure generation failed — continuing"
fi

log "Pipeline complete in $(( ($(date +%s) - started) / 60 )) min"
