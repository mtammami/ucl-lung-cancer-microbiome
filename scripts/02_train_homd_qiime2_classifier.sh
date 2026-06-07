#!/usr/bin/env bash



# Train a QIIME 2 Naive Bayes classifier using HOMD/eHOMD 16S reference sequences.


HOMD_FASTA="${1:-databases/homd/homd_refseq.fasta}"
HOMD_TAXONOMY="${2:-databases/homd/homd_taxonomy_qiime.txt}"
OUTPUT_DIR="${3:-databases/homd/classifier}"
MODE="${4:-v3v4}"
THREADS="${5:-8}"

FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
REVERSE_PRIMER="GACTACHVGGGTATCTAATCC"

MIN_LENGTH=250
MAX_LENGTH=550

mkdir -p "$OUTPUT_DIR"

CLEAN_FASTA="${OUTPUT_DIR}/homd_refseq_cleaned.fasta"

REF_SEQS_QZA="${OUTPUT_DIR}/homd_refseq.qza"
REF_TAX_QZA="${OUTPUT_DIR}/homd_taxonomy.qza"

EXTRACTED_READS_QZA="${OUTPUT_DIR}/homd_${MODE}_reads.qza"
CLASSIFIER_QZA="${OUTPUT_DIR}/homd_${MODE}_classifier.qza"

PREDICTED_TAX_QZA="${OUTPUT_DIR}/homd_${MODE}_training_predictions.qza"
PREDICTED_TAX_QZV="${OUTPUT_DIR}/homd_${MODE}_training_predictions.qzv"

if ! command -v qiime >/dev/null 2>&1; then
  echo "[ERROR] qiime command not found. Activate your QIIME 2 environment first."
  exit 1
fi

if [[ ! -f "$HOMD_FASTA" ]]; then
  echo "[ERROR] HOMD FASTA file not found:"
  echo "$HOMD_FASTA"
  exit 1
fi

if [[ ! -f "$HOMD_TAXONOMY" ]]; then
  echo "[ERROR] HOMD taxonomy file not found:"
  echo "$HOMD_TAXONOMY"
  exit 1
fi

if [[ "$MODE" != "v3v4" && "$MODE" != "full_length" ]]; then
  echo "[ERROR] MODE must be either: v3v4 or full_length"
  exit 1
fi

echo "========================================================="
echo "Training HOMD QIIME 2 classifier"
echo "========================================================="
echo "HOMD FASTA:       $HOMD_FASTA"
echo "HOMD taxonomy:    $HOMD_TAXONOMY"
echo "Output directory: $OUTPUT_DIR"
echo "Mode:             $MODE"
echo "Threads:          $THREADS"
echo "Forward primer:   $FORWARD_PRIMER"
echo "Reverse primer:   $REVERSE_PRIMER"
echo "========================================================="

echo "[STEP 1] Cleaning HOMD FASTA file..."

awk '
  /^>/ {
    print
    next
  }
  {
    gsub(/[-.[:space:]]/, "", $0)
    gsub(/u/, "t", $0)
    gsub(/U/, "T", $0)
    print toupper($0)
  }
' "$HOMD_FASTA" > "$CLEAN_FASTA"

echo "[DONE] Cleaned FASTA:"
echo "$CLEAN_FASTA"

echo "[STEP 2] Importing HOMD reference sequences..."

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "$CLEAN_FASTA" \
  --output-path "$REF_SEQS_QZA"

echo "[DONE] Reference sequences artifact:"
echo "$REF_SEQS_QZA"

echo "[STEP 3] Importing HOMD taxonomy..."

FIRST_LINE="$(head -n 1 "$HOMD_TAXONOMY")"

if echo "$FIRST_LINE" | grep -Eiq 'Feature|Taxon|taxonomy'; then
  TAX_FORMAT="TSVTaxonomyFormat"
else
  TAX_FORMAT="HeaderlessTSVTaxonomyFormat"
fi

echo "[INFO] Detected taxonomy format: $TAX_FORMAT"

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format "$TAX_FORMAT" \
  --input-path "$HOMD_TAXONOMY" \
  --output-path "$REF_TAX_QZA"

echo "[DONE] Taxonomy artifact:"
echo "$REF_TAX_QZA"

if [[ "$MODE" == "v3v4" ]]; then
  echo "[STEP 4] Extracting V3-V4 reads from HOMD reference sequences..."

  qiime feature-classifier extract-reads \
    --i-sequences "$REF_SEQS_QZA" \
    --p-f-primer "$FORWARD_PRIMER" \
    --p-r-primer "$REVERSE_PRIMER" \
    --p-min-length "$MIN_LENGTH" \
    --p-max-length "$MAX_LENGTH" \
    --p-n-jobs "$THREADS" \
    --o-reads "$EXTRACTED_READS_QZA"

  TRAINING_READS_QZA="$EXTRACTED_READS_QZA"

  echo "[DONE] Extracted V3-V4 reads:"
  echo "$EXTRACTED_READS_QZA"

else
  echo "[STEP 4] Skipping primer extraction."
  echo "[INFO] Training classifier with full-length HOMD sequences."

  TRAINING_READS_QZA="$REF_SEQS_QZA"
fi

echo "[STEP 5] Training QIIME 2 Naive Bayes classifier..."

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "$TRAINING_READS_QZA" \
  --i-reference-taxonomy "$REF_TAX_QZA" \
  --o-classifier "$CLASSIFIER_QZA"

echo "[DONE] HOMD classifier:"
echo "$CLASSIFIER_QZA"

echo "[STEP 6] Testing classifier on HOMD training reads..."

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_QZA" \
  --i-reads "$TRAINING_READS_QZA" \
  --o-classification "$PREDICTED_TAX_QZA" \
  --p-n-jobs "$THREADS" \
  --verbose

qiime metadata tabulate \
  --m-input-file "$PREDICTED_TAX_QZA" \
  --o-visualization "$PREDICTED_TAX_QZV"

echo "[DONE] Training-set classification:"
echo "$PREDICTED_TAX_QZA"

echo "[DONE] Training-set classification visualization:"
echo "$PREDICTED_TAX_QZV"

echo "========================================================="
echo "HOMD classifier training completed."
echo "========================================================="
echo "Main classifier:"
echo "  $CLASSIFIER_QZA"
echo ""
echo "Reference artifacts:"
echo "  $REF_SEQS_QZA"
echo "  $REF_TAX_QZA"
echo ""
echo "Validation output:"
echo "  $PREDICTED_TAX_QZV"
echo "========================================================="