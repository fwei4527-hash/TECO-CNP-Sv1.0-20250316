#!/usr/bin/env bash
set -Eeuo pipefail

# ============================================================
# TECO-CNP three-stage simulation
#
# Stage 1: 1934-1983 forcing repeated for equilibrium spin-up
# Stage 2: 1984-2021 forest establishment and historical simulation
# Stage 3: 2022-2024 P0/P25/P50/P100 experiments
#
# Initialstate_heshan.csv:
# row 1 = nine carbon pools, g C m-2
# row 2 = nine C:N ratios
# row 3 = nine C:P ratios
# ============================================================

MAIN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${MAIN_DIR}/.." && pwd)"

INPUT_DIR="${ROOT_DIR}/input"
OUTPUT_DIR="${ROOT_DIR}/output"
FORCING_DIR="${INPUT_DIR}/three_stage_forcing"

ACTIVE_FORCING="${INPUT_DIR}/TECO forcing 2022-2024.txt"

# Fortran实际读取的文件
ACTIVE_INITIAL="${INPUT_DIR}/Initialstate_heshan.csv"

# 阶段1使用的低碳库初始文件
STAGE1_SEED_INITIAL="${INPUT_DIR}/Initialstate_heshan zero.csv"

STAGE1_FORCING="${FORCING_DIR}/TECO_阶段1_预热_1934_1983.txt"
STAGE2_FORCING="${FORCING_DIR}/TECO_阶段2_历史_1984_2021.txt"
STAGE3_FORCING="${FORCING_DIR}/TECO_阶段3_公共气象_2022_2024.txt"

# 50-year forcing × 100 cycles = 5000 simulated years.
# Do not set this to 5000 unless you intend to simulate 250000 years.
N_SPINUP="${N_SPINUP:-100}"

SPIN_REL_TOL="${SPIN_REL_TOL:-1e-4}"
SPIN_ABS_TOL="${SPIN_ABS_TOL:-1e-4}"

# 0.1 g C m-2 = 0.001 t C ha-1, consistent with the small
# initial plant pools used in the GDAY historical establishment run.
PLANT_INITIAL_C="${PLANT_INITIAL_C:-0.1}"

TIME_TAG="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUTPUT_DIR}/three_stage_${TIME_TAG}"
BACKUP_DIR="${RUN_DIR}/input_backup"
BUILD_DIR="${MAIN_DIR}/build_three_stage"

SPIN_MODEL_OUTPUT="${OUTPUT_DIR}/NDSPINUP/teco_cnp"
SIM_MODEL_OUTPUT="${OUTPUT_DIR}/sim/teco_cnp"

mkdir -p "${RUN_DIR}" "${BACKUP_DIR}" "${BUILD_DIR}"

require_file() {
    if [[ ! -s "$1" ]]; then
        echo "ERROR: missing or empty file: $1"
        exit 1
    fi
}

require_file "${ACTIVE_INITIAL}"
require_file "${ACTIVE_FORCING}"
require_file "${STAGE1_FORCING}"
require_file "${STAGE2_FORCING}"
require_file "${STAGE3_FORCING}"
require_file "${STAGE1_SEED_INITIAL}"

# Preserve the user's current input files.
cp -p "${ACTIVE_INITIAL}" \
    "${BACKUP_DIR}/Initialstate_heshan.csv"

cp -p "${ACTIVE_FORCING}" \
    "${BACKUP_DIR}/TECO_forcing_original.txt"

restore_inputs() {
    cp -f "${BACKUP_DIR}/Initialstate_heshan.csv" \
        "${ACTIVE_INITIAL}"

    cp -f "${BACKUP_DIR}/TECO_forcing_original.txt" \
        "${ACTIVE_FORCING}"
}

trap restore_inputs EXIT

# The initial file present when the script starts is used only as
# the numerical seed for Stage 1.
STAGE1_INITIAL="${RUN_DIR}/Initialstate_stage1_seed.csv"
cp -p "${STAGE1_SEED_INITIAL}" "${STAGE1_INITIAL}"

# ------------------------------------------------------------
# Compile once
# ------------------------------------------------------------
echo "Compiling TECO-CNP..."

cd "${MAIN_DIR}"

SOURCES=(
    FileSize.f90
    ParasModule.f90
    SASpinUp.f90
    LIMITATION.f90
    NPUptakeDemand.f90
    NPDynamic.f90
    MCMC.f90
    TECO_CNP_main.f90
)

OBJECTS=()

for source_file in "${SOURCES[@]}"; do
    object_file="${BUILD_DIR}/${source_file%.f90}.o"

    gfortran -O2 -ffree-form \
        -J"${BUILD_DIR}" \
        -I"${BUILD_DIR}" \
        -c "${source_file}" \
        -o "${object_file}"

    OBJECTS+=("${object_file}")
done

gfortran -o "${MAIN_DIR}/teco_cnp.exe" \
    "${OBJECTS[@]}" \
    -llapacke -llapack -lblas

require_file "${MAIN_DIR}/teco_cnp.exe"

# ------------------------------------------------------------
# Output handling
# ------------------------------------------------------------
prepare_output_directory() {
    local model_output="$1"
    local backup_name="$2"

    if [[ -d "${model_output}" ]]; then
        mv "${model_output}" \
            "${RUN_DIR}/${backup_name}_preexisting"
    fi

    mkdir -p "${model_output}"
}

save_model_output() {
    local model_output="$1"
    local destination="$2"

    require_file "${model_output}/teco_8cpools_cnp.csv"
    require_file "${model_output}/CNP_Ncycle_simu_ratios_cnp.csv"
    require_file "${model_output}/CNP_Pcycle_simu_ratios_cnp.csv"

    mv "${model_output}" "${destination}"
}

# ------------------------------------------------------------
# Construct the next Initialstate file from model output.
#
# reset_plant = 1:
#   reset QC(1:4), retain spun-up soil QC(5:9).
#
# reset_plant = 0:
#   retain all nine carbon pools.
# ------------------------------------------------------------
make_initial_state() {
    local result_dir="$1"
    local output_file="$2"
    local reset_plant="$3"

    python3 - \
        "${result_dir}" \
        "${output_file}" \
        "${reset_plant}" \
        "${PLANT_INITIAL_C}" <<'PY'
import csv
import math
import pathlib
import sys

result_dir = pathlib.Path(sys.argv[1])
output_file = pathlib.Path(sys.argv[2])
reset_plant = int(sys.argv[3])
plant_initial = float(sys.argv[4])

files = {
    "carbon": result_dir / "teco_8cpools_cnp.csv",
    "cn": result_dir / "CNP_Ncycle_simu_ratios_cnp.csv",
    "cp": result_dir / "CNP_Pcycle_simu_ratios_cnp.csv",
}

def read_last_row(path):
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"Missing output file: {path}")

    last = None
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            values = [item.strip() for item in row if item.strip()]
            if values:
                last = [float(item) for item in values]

    if last is None or len(last) < 9:
        raise SystemExit(f"Invalid final row in {path}")

    last = last[:9]

    if not all(math.isfinite(value) for value in last):
        raise SystemExit(f"Non-finite state found in {path}")

    return last

carbon = read_last_row(files["carbon"])
cn_ratio = read_last_row(files["cn"])
cp_ratio = read_last_row(files["cp"])

if reset_plant:
    carbon[0:4] = [plant_initial] * 4

if any(value <= 0 for value in carbon):
    raise SystemExit("All carbon pools must be greater than zero.")

if any(value <= 0 for value in cn_ratio):
    raise SystemExit("All C:N ratios must be greater than zero.")

if any(value <= 0 for value in cp_ratio):
    raise SystemExit("All C:P ratios must be greater than zero.")

output_file.parent.mkdir(parents=True, exist_ok=True)

with output_file.open("w", encoding="utf-8", newline="\n") as handle:
    for values in (carbon, cn_ratio, cp_ratio):
        handle.write(",".join(f"{value:.10g}" for value in values))
        handle.write("\n")

print(f"Created initial state: {output_file}")
PY
}

# ------------------------------------------------------------
# Check the final two spin-up cycles.
#
# Carbon: columns 5-13 in spinup_LoopVariables
# Nitrogen: first nine columns in spinup_LoopVariablesNNN
# Phosphorus: first nine columns in spinup_LoopVariablesPPP
# ------------------------------------------------------------
check_spinup_convergence() {
    local result_dir="$1"

    python3 - \
        "${result_dir}" \
        "${SPIN_REL_TOL}" \
        "${SPIN_ABS_TOL}" <<'PY'
import csv
import math
import pathlib
import sys

result_dir = pathlib.Path(sys.argv[1])
relative_tolerance = float(sys.argv[2])
absolute_tolerance = float(sys.argv[3])

files = {
    "C": (
        result_dir / "spinup_LoopVariables_cnp.csv",
        slice(4, 13)
    ),
    "N": (
        result_dir / "spinup_LoopVariablesNNN_cnp.csv",
        slice(0, 9)
    ),
    "P": (
        result_dir / "spinup_LoopVariablesPPP_cnp.csv",
        slice(0, 9)
    ),
}

failed = []

for element, (path, selected_columns) in files.items():
    if not path.exists():
        raise SystemExit(f"Missing spin-up check file: {path}")

    rows = []

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            values = [item.strip() for item in row if item.strip()]
            if values:
                rows.append([float(item) for item in values])

    if len(rows) < 2:
        raise SystemExit(f"Not enough spin-up cycles in {path}")

    previous = rows[-2][selected_columns]
    current = rows[-1][selected_columns]

    for pool_number, (old, new) in enumerate(
        zip(previous, current), start=1
    ):
        scale = max(abs(old), abs(new), 1.0)
        threshold = absolute_tolerance + relative_tolerance * scale
        difference = abs(new - old)

        if (
            not math.isfinite(old)
            or not math.isfinite(new)
            or difference > threshold
        ):
            failed.append(
                (
                    element,
                    pool_number,
                    old,
                    new,
                    difference,
                    threshold,
                )
            )

if failed:
    print("Spin-up has not converged:")
    for item in failed[:30]:
        element, pool, old, new, difference, threshold = item
        print(
            f"{element} pool {pool}: "
            f"previous={old:.8g}, current={new:.8g}, "
            f"difference={difference:.4g}, "
            f"threshold={threshold:.4g}"
        )
    raise SystemExit(
        "Increase N_SPINUP or inspect unstable pools."
    )

print("Spin-up convergence check passed.")
PY
}

# ============================================================
# Stage 1: pre-1984 equilibrium spin-up
# ============================================================
echo "Stage 1: 1934-1983 equilibrium spin-up"

cp -f "${STAGE1_FORCING}" "${ACTIVE_FORCING}"
cp -f "${STAGE1_INITIAL}" "${ACTIVE_INITIAL}"

prepare_output_directory \
    "${SPIN_MODEL_OUTPUT}" \
    "stage1_spinup"

"${MAIN_DIR}/teco_cnp.exe" \
    1934 1983 \
    3 \
    0 \
    1 \
    "${N_SPINUP}" \
    0 \
    0 \
    0 \
    > "${RUN_DIR}/stage1_spinup.log" 2>&1

save_model_output \
    "${SPIN_MODEL_OUTPUT}" \
    "${RUN_DIR}/stage1_spinup"

check_spinup_convergence "${RUN_DIR}/stage1_spinup"

# Build the 1984 initial state:
# plant pools small, soil pools inherited from spin-up.
STAGE2_INITIAL="${RUN_DIR}/Initialstate_stage2_1984.csv"

make_initial_state \
    "${RUN_DIR}/stage1_spinup" \
    "${STAGE2_INITIAL}" \
    1

# ============================================================
# Stage 2: 1984-2021 historical simulation
# ============================================================
echo "Stage 2: 1984-2021 historical simulation"

cp -f "${STAGE2_FORCING}" "${ACTIVE_FORCING}"
cp -f "${STAGE2_INITIAL}" "${ACTIVE_INITIAL}"

prepare_output_directory \
    "${SIM_MODEL_OUTPUT}" \
    "stage2_historical"

"${MAIN_DIR}/teco_cnp.exe" \
    1984 2021 \
    3 \
    0 \
    0 \
    1 \
    0 \
    0 \
    0 \
    > "${RUN_DIR}/stage2_historical.log" 2>&1

save_model_output \
    "${SIM_MODEL_OUTPUT}" \
    "${RUN_DIR}/stage2_historical"

# Build the common state at the beginning of 2022.
STAGE3_INITIAL="${RUN_DIR}/Initialstate_stage3_2022.csv"

make_initial_state \
    "${RUN_DIR}/stage2_historical" \
    "${STAGE3_INITIAL}" \
    0

# ============================================================
# Stage 3: 2022-2024 P-addition treatments
# ============================================================
echo "Stage 3: 2022-2024 P-addition treatments"

cp -f "${STAGE3_FORCING}" "${ACTIVE_FORCING}"

for P_RATE in 0 25 50 100; do

    echo "Running P${P_RATE}"

    cp -f "${STAGE3_INITIAL}" "${ACTIVE_INITIAL}"

    prepare_output_directory \
        "${SIM_MODEL_OUTPUT}" \
        "stage3_P${P_RATE}"

    if [[ "${P_RATE}" -eq 0 ]]; then
        NP_ADDITION=0
    else
        NP_ADDITION=1
    fi

    "${MAIN_DIR}/teco_cnp.exe" \
        2022 2024 \
        3 \
        0 \
        0 \
        1 \
        0 \
        "${NP_ADDITION}" \
        "${P_RATE}" \
        > "${RUN_DIR}/stage3_P${P_RATE}.log" 2>&1

    save_model_output \
        "${SIM_MODEL_OUTPUT}" \
        "${RUN_DIR}/stage3_P${P_RATE}"

done

echo "All three stages completed."
echo "Results: ${RUN_DIR}"