#!/usr/bin/env bash
set -euo pipefail

# Low-memory split-and-filter workflow for filtering NCBI NR with seqkit.
# Configure paths via environment variables or edit the defaults below.
NR_GZ="${NR_GZ:-nr.gz}"
SEQIDS_PREFIX="${SEQIDS_PREFIX:-keep_seqids_part_}"
OUT_DIR="${OUT_DIR:-filtered_nr}"
SEQKIT_BIN="${SEQKIT_BIN:-seqkit}"
THREADS="${THREADS:-16}"
PIGZ_THREADS="${PIGZ_THREADS:-8}"
CLEANUP_PARTS="${CLEANUP_PARTS:-0}"

for cmd in pigz "${SEQKIT_BIN}"; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: Required command not found on PATH: ${cmd}" >&2
    exit 1
  fi
done

if [[ ! -f "${NR_GZ}" ]]; then
  echo "ERROR: Input NR archive not found: ${NR_GZ}" >&2
  echo "Set NR_GZ=/path/to/nr.gz before running this script." >&2
  exit 1
fi

shopt -s nullglob
part_files=("${SEQIDS_PREFIX}"*)
shopt -u nullglob

if [[ "${#part_files[@]}" -eq 0 ]]; then
  echo "ERROR: No chunk files matched '${SEQIDS_PREFIX}*'." >&2
  echo "Split your accession list first, e.g.: split -n l/3 -d keep_seqids.txt ${SEQIDS_PREFIX}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"
LOG_FILE="${OUT_DIR}/filter_log.txt"

echo "Process started at $(date)" | tee "${LOG_FILE}"
echo "Input NR: ${NR_GZ}" | tee -a "${LOG_FILE}"
echo "Chunks found: ${#part_files[@]}" | tee -a "${LOG_FILE}"

filtered_parts=()
pass=0
for part_file in "${part_files[@]}"; do
  pass=$((pass + 1))
  part_name="$(basename "${part_file}")"
  out_part="${OUT_DIR}/nr_filtered_${part_name}.faa.gz"

  echo "Starting Pass ${pass}/${#part_files[@]} (${part_name}) at $(date)..." | tee -a "${LOG_FILE}"
  pigz -dc "${NR_GZ}" | \
    "${SEQKIT_BIN}" grep -f "${part_file}" --id-regexp '^([^\s]+)' -j "${THREADS}" | \
    pigz -p "${PIGZ_THREADS}" > "${out_part}"

  filtered_parts+=("${out_part}")
done

MERGED_OUT="${OUT_DIR}/nr_filtered.faa.gz"
echo "Merging ${#filtered_parts[@]} filtered parts into ${MERGED_OUT} at $(date)..." | tee -a "${LOG_FILE}"
cat "${filtered_parts[@]}" > "${MERGED_OUT}"

if [[ "${CLEANUP_PARTS}" == "1" ]]; then
  rm -f "${filtered_parts[@]}"
fi

echo "Completed successfully at $(date)!" | tee -a "${LOG_FILE}"
