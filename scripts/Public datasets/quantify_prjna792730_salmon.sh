#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${MPP_PUBLIC_ROOT:-inputs/public}"
FASTQ_DIR="$BASE_DIR/raw/PRJNA792730/fastq"
REF_DIR="$BASE_DIR/reference/PRJNA792730/gencode_v44"
INDEX_DIR="$REF_DIR/salmon_index"
QUANT_DIR="$BASE_DIR/quant/PRJNA792730"
RESULT_DIR="$BASE_DIR/results/PRJNA792730"
MANIFEST="$BASE_DIR/results/PRJNA792730_ena_read_run_with_md5.tsv"
SALMON="${SALMON:-salmon}"
TRANSCRIPTS="$REF_DIR/gencode.v44.transcripts.fa.gz"

THREADS="${THREADS:-8}"

declare -a RUNS=(
  SRR17345233 SRR17345232 SRR17345231
  SRR17345230 SRR17345229 SRR17345228
  SRR17345227 SRR17345226 SRR17345225
)
declare -a SAMPLES=(
  C1-0h C2-0h C3-0h
  T1-12h T2-12h T3-12h
  T1-24h T2-24h T3-24h
)

mkdir -p "$INDEX_DIR" "$QUANT_DIR" "$RESULT_DIR"

if [[ ! -x "$SALMON" ]]; then
  printf 'Salmon executable not found: %s\n' "$SALMON" >&2
  exit 1
fi
for required in "$MANIFEST" "$TRANSCRIPTS"; do
  if [[ ! -s "$required" ]]; then
    printf 'Required input is missing or empty: %s\n' "$required" >&2
    exit 1
  fi
done
gzip -t "$TRANSCRIPTS"

printf 'run\tsample\tread\texpected_bytes\tobserved_bytes\texpected_md5\tobserved_md5\tgzip_test\n' \
  > "$RESULT_DIR/PRJNA792730_fastq_prequant_qc.tsv"

for idx in "${!RUNS[@]}"; do
  run="${RUNS[$idx]}"
  sample="${SAMPLES[$idx]}"
  manifest_row="$(awk -F '\t' -v run="$run" 'NR > 1 && $1 == run {print; exit}' "$MANIFEST")"
  if [[ -z "$manifest_row" ]]; then
    printf 'Run is absent from the ENA manifest: %s\n' "$run" >&2
    exit 1
  fi
  fastq_bytes="$(printf '%s\n' "$manifest_row" | awk -F '\t' '{print $4}')"
  fastq_md5="$(printf '%s\n' "$manifest_row" | awk -F '\t' '{print $5}')"
  IFS=';' read -r expected_r1 expected_r2 <<< "$fastq_bytes"
  IFS=';' read -r expected_md5_r1 expected_md5_r2 <<< "$fastq_md5"
  for read in 1 2; do
    file="$FASTQ_DIR/${run}_${read}.fastq.gz"
    if [[ ! -f "$file" ]]; then
      printf 'FASTQ is missing: %s\n' "$file" >&2
      exit 1
    fi
    expected_var="expected_r${read}"
    expected="${!expected_var}"
    observed="$(stat -f '%z' "$file")"
    if [[ "$observed" != "$expected" ]]; then
      printf 'FASTQ size mismatch: %s observed=%s expected=%s\n' \
        "$file" "$observed" "$expected" >&2
      exit 1
    fi
    gzip -t "$file"
    expected_md5_var="expected_md5_r${read}"
    expected_md5="${!expected_md5_var}"
    observed_md5="$(md5 -q "$file")"
    if [[ "$observed_md5" != "$expected_md5" ]]; then
      printf 'FASTQ MD5 mismatch: %s observed=%s expected=%s\n' \
        "$file" "$observed_md5" "$expected_md5" >&2
      exit 1
    fi
    printf '%s\t%s\tR%s\t%s\t%s\t%s\t%s\tPASS\n' \
      "$run" "$sample" "$read" "$expected" "$observed" "$expected_md5" "$observed_md5" \
      >> "$RESULT_DIR/PRJNA792730_fastq_prequant_qc.tsv"
  done
done

if [[ ! -s "$INDEX_DIR/versionInfo.json" ]]; then
  if find "$INDEX_DIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    mv "$INDEX_DIR" "$INDEX_DIR.incomplete.$(date +%s)"
    mkdir -p "$INDEX_DIR"
  fi
  "$SALMON" index \
    --transcripts "$TRANSCRIPTS" \
    --index "$INDEX_DIR" \
    --gencode \
    --threads "$THREADS"
fi

for idx in "${!RUNS[@]}"; do
  run="${RUNS[$idx]}"
  sample="${SAMPLES[$idx]}"
  sample_dir="$QUANT_DIR/$sample"
  if [[ -s "$sample_dir/quant.sf" && -s "$sample_dir/aux_info/meta_info.json" ]]; then
    printf 'Quantification already complete: %s\n' "$sample"
    continue
  fi
  if [[ -e "$sample_dir.incomplete" ]]; then
    mv "$sample_dir.incomplete" "$sample_dir.incomplete.$(date +%s)"
  fi
  "$SALMON" quant \
    --index "$INDEX_DIR" \
    --libType A \
    --mates1 "$FASTQ_DIR/${run}_1.fastq.gz" \
    --mates2 "$FASTQ_DIR/${run}_2.fastq.gz" \
    --validateMappings \
    --seqBias \
    --gcBias \
    --threads "$THREADS" \
    --output "$sample_dir.incomplete"
  mv "$sample_dir.incomplete" "$sample_dir"
done

Rscript "$BASE_DIR/build_prjna792730_gene_counts.R" \
  "$REF_DIR/gencode.v44.annotation.gtf.gz" \
  "$QUANT_DIR" \
  "$RESULT_DIR/PRJNA792730_salmon_gene_counts.tsv" \
  "$RESULT_DIR/PRJNA792730_salmon_qc.tsv"

Rscript "$BASE_DIR/analyze_prjna792730_counts.R" \
  "$RESULT_DIR/PRJNA792730_salmon_gene_counts.tsv" \
  "$RESULT_DIR/analysis"

printf 'PRJNA792730 quantification and prespecified analysis completed.\n'
