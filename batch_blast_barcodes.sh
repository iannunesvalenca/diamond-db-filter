#!/usr/bin/env bash
set -euo pipefail

# Batch DIAMOND BLASTx for barcoded FASTQ files with NCBI taxonomy summaries.
# Configure via positional arguments or environment variables:
#   Usage: bash batch_blast_barcodes.sh [INPUT_DIR] [DB_PATH] [OUTPUT_DIR]
#   Env:   TAXDMP_DIR, THREADS, EVALUE, MAX_TARGET_SEQS, FRAMESHIFT

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_DIR="${1:-${INPUT_DIR:-}}"
DB="${2:-${DB:-nr_filtered.dmnd}}"
OUTPUT_DIR="${3:-${OUTPUT_DIR:-diamond_reports_viral}}"
TAXDMP_DIR="${TAXDMP_DIR:-taxdmp}"
THREADS="${THREADS:-16}"
EVALUE="${EVALUE:-0.001}"
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-1}"
FRAMESHIFT="${FRAMESHIFT:-15}"

if [[ -z "${INPUT_DIR}" || ! -d "${INPUT_DIR}" ]]; then
  echo "ERROR: Input FASTQ directory not found: '${INPUT_DIR}'" >&2
  echo "Usage: bash $0 <INPUT_DIR> [DB_PATH] [OUTPUT_DIR]" >&2
  exit 1
fi

if [[ ! -f "${DB}" ]]; then
  echo "ERROR: DIAMOND database not found: ${DB}" >&2
  exit 1
fi

for cmd in diamond python3; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: Required command not found on PATH: ${cmd}" >&2
    exit 1
  fi
done

mkdir -p "${OUTPUT_DIR}"

echo "Starting batch DIAMOND BLASTx..."
echo "  Input dir:  ${INPUT_DIR}"
echo "  Database:   ${DB}"
echo "  Output dir: ${OUTPUT_DIR}"
echo "  Params:     --very-sensitive --evalue ${EVALUE} --max-target-seqs ${MAX_TARGET_SEQS} --frameshift ${FRAMESHIFT} --threads ${THREADS}"

shopt -s nullglob
fastq_files=("${INPUT_DIR}"/*.fastq "${INPUT_DIR}"/*.fastq.gz)
shopt -u nullglob

if [[ "${#fastq_files[@]}" -eq 0 ]]; then
  echo "No .fastq or .fastq.gz files found in ${INPUT_DIR}."
  exit 0
fi

count=0
for fastq in "${fastq_files[@]}"; do
  base="$(basename "${fastq}")"
  filename="${base%.gz}"
  filename="${filename%.fastq}"
  output_file="${OUTPUT_DIR}/${filename}_viral_report.tsv"

  echo "---------------------------------------------------"
  echo "Processing: ${filename}"

  if diamond blastx \
      -q "${fastq}" \
      -d "${DB}" \
      -o "${output_file}" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms \
      --evalue "${EVALUE}" \
      --max-target-seqs "${MAX_TARGET_SEQS}" \
      --very-sensitive \
      --frameshift "${FRAMESHIFT}" \
      --threads "${THREADS}"; then
    echo "BLASTx complete. Running taxonomy summarizer..."
    python3 "${SCRIPT_DIR}/summarize_diamond.py" \
      "${output_file}" \
      --nodes "${TAXDMP_DIR}/nodes.dmp" \
      --names "${TAXDMP_DIR}/names.dmp" \
      --merged "${TAXDMP_DIR}/merged.dmp"
  else
    echo "ERROR: DIAMOND BLASTx failed for ${filename}" >&2
  fi
  count=$((count + 1))
done

echo "Processing complete (${count} files)."
