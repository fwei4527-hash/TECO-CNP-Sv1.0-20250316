#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL_ROOT="${1:-$HOME/fangwei/FangxiuWan-TECO-CNP-Sv1.0-edff198/TECO-CNP Sv1.0 20250316}"
TRAJECTORIES="${TRAJECTORIES:-20}"
WORKERS="${WORKERS:-2}"
SEED="${SEED:-20260921}"
RESULT_DIR="$MODEL_ROOT/output/sensitivity/TECO_C_P_coupling_Morris"

python3 "$SCRIPT_DIR/run_teco_morris.py" \
    --model-root "$MODEL_ROOT" \
    --executable "$MODEL_ROOT/main/teco_sensitivity.exe" \
    --candidates "$SCRIPT_DIR/teco_sensitivity_candidates.csv" \
    --metric-groups "$SCRIPT_DIR/metric_groups.json" \
    --results "$RESULT_DIR" \
    --start-year 2021 \
    --end-year 2024 \
    --trajectories "$TRAJECTORIES" \
    --levels 6 \
    --seed "$SEED" \
    --workers "$WORKERS" \
    --top-per-process 3
