#!/usr/bin/env bash
set -Eeuo pipefail

# Launch five LocalColabFold workers, one per GPU.
#
# Usage:
#   bash pipeline/run_colabfold_cluster.sh a3m_batches colabfold_outputs logs_colabfold

A3M_ROOT="${1:-a3m_batches}"
OUT_ROOT="${2:-colabfold_outputs}"
LOG_ROOT="${3:-logs_colabfold}"

NUM_GPU=5
COLABFOLD_BIN="${COLABFOLD_BIN:-colabfold_batch}"
NUM_RECYCLE=3
NUM_RELAX=0
NUM_MODELS=1
MODEL_ORDER=1
STOP_AT_SCORE=100

export XLA_PYTHON_CLIENT_MEM_FRACTION="${XLA_PYTHON_CLIENT_MEM_FRACTION:-0.9}"
export TF_FORCE_UNIFIED_MEMORY="${TF_FORCE_UNIFIED_MEMORY:-1}"
export XLA_PYTHON_CLIENT_PREALLOCATE="${XLA_PYTHON_CLIENT_PREALLOCATE:-true}"

mkdir -p "${OUT_ROOT}" "${LOG_ROOT}"

if ! command -v "${COLABFOLD_BIN}" >/dev/null 2>&1; then
  echo "[ERROR] colabfold_batch not found in PATH: ${COLABFOLD_BIN}" >&2
  exit 127
fi

for gpu in $(seq 0 $((NUM_GPU - 1))); do
  batch_dir="${A3M_ROOT}/batch_${gpu}"
  if [[ ! -d "${batch_dir}" ]]; then
    echo "[ERROR] Missing batch directory: ${batch_dir}" >&2
    exit 2
  fi
  count="$(find "${batch_dir}" -maxdepth 1 -type f -name '*.a3m' | wc -l | tr -d ' ')"
  if [[ "${count}" == "0" ]]; then
    echo "[ERROR] No .a3m files found in ${batch_dir}" >&2
    exit 2
  fi
  echo "[INFO] batch_${gpu}: ${count} A3M files"
done

declare -a PIDS=()

cleanup() {
  local exit_code=$?
  if [[ "${exit_code}" -ne 0 ]]; then
    echo "[ERROR] Terminating ColabFold workers" >&2
    for pid in "${PIDS[@]:-}"; do
      if kill -0 "${pid}" >/dev/null 2>&1; then
        kill "${pid}" >/dev/null 2>&1 || true
      fi
    done
  fi
  exit "${exit_code}"
}
trap cleanup EXIT INT TERM

for gpu in $(seq 0 $((NUM_GPU - 1))); do
  batch_dir="${A3M_ROOT}/batch_${gpu}"
  out_dir="${OUT_ROOT}/batch_${gpu}"
  log_file="${LOG_ROOT}/batch_${gpu}.log"
  mkdir -p "${out_dir}"

  echo "[INFO] Starting GPU ${gpu}: ${batch_dir} -> ${out_dir}"
  (
    set -Eeuo pipefail
    export CUDA_VISIBLE_DEVICES="${gpu}"
    export JAX_PLATFORM_NAME="gpu"
    "${COLABFOLD_BIN}" \
      --msa-mode custom \
      --num-recycle "${NUM_RECYCLE}" \
      --num-relax "${NUM_RELAX}" \
      --num-models "${NUM_MODELS}" \
      --model-order "${MODEL_ORDER}" \
      --stop-at-score "${STOP_AT_SCORE}" \
      "${batch_dir}" \
      "${out_dir}"
  ) >"${log_file}" 2>&1 &
  PIDS+=("$!")
done

failed=0
for idx in "${!PIDS[@]}"; do
  if wait "${PIDS[$idx]}"; then
    echo "[INFO] Worker ${idx} finished"
  else
    echo "[ERROR] Worker ${idx} failed. See ${LOG_ROOT}/batch_${idx}.log" >&2
    failed=1
  fi
done

trap - EXIT INT TERM

if [[ "${failed}" -ne 0 ]]; then
  exit 1
fi

echo "[INFO] All ColabFold workers completed"
