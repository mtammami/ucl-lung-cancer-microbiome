#!/usr/bin/env bash


# QIIME 2 paired-end 16S V3-V4 pipeline
#
# Steps:
#   1. Import paired-end FASTQ files
#   2. Trim V3-V4 primers using cutadapt
#   3. Summarize demultiplexed reads
#   4. DADA2 denoising
#   5. DADA2 statistics summary
#   6. Feature table and representative sequence summaries
#   7. Taxonomic classification using HOMD classifier
#   8. Phylogenetic tree construction
#
# Usage:
#   bash scripts/09_qiime2_miseq_v3v4_pipeline.sh \
#     <input_fastq_directory> \
#     <output_directory> \
#     <HOMD.qza> \
#     <threads>
#
# Example:
#   bash scripts/09_qiime2_miseq_v3v4_pipeline.sh \
#     data/raw/miseq \
#     results/qiime2 \
#     databases/HOMD-classifier.qza \
#     20

# ---------------- User arguments ----------------

INPUT_DIR="${1:-data/raw/miseq}"
OUTPUT_DIR="${2:-results/qiime2}"
CLASSIFIER="${3:-databases/HOMD-classifier.qza}"
THREADS="${4:-20}"

# ---------------- Parameters ----------------

# Correct 16S V3-V4 primers
FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
REVERSE_PRIMER="GACTACHVGGGTATCTAATCC"

TRIM_LEFT_F=0
TRIM_LEFT_R=0
TRUNC_LEN_F=240
TRUNC_LEN_R=200

# ---------------- Output files ----------------

DEMUX_QZA="${OUTPUT_DIR}/demux-paired-end.qza"
TRIMMED_QZA="${OUTPUT_DIR}/primer-trimmed-paired-end.qza"
TRIMMED_QZV="${OUTPUT_DIR}/primer-trimmed-paired-end.qzv"

DADA2_DIR="${OUTPUT_DIR}/dada2_denoising_${TRUNC_LEN_F}_${TRUNC_LEN_R}"
DADA2_TABLE="${DADA2_DIR}/table.qza"
DADA2_REP_SEQS="${DADA2_DIR}/representative_sequences.qza"
DADA2_STATS="${DADA2_DIR}/denoising_stats.qza"

DADA2_STATS_QZV="${OUTPUT_DIR}/denoising-stats.qzv"
FEATURE_TABLE_QZV="${OUTPUT_DIR}/feature-table-summary.qzv"
REP_SEQS_QZV="${OUTPUT_DIR}/rep-seqs.qzv"

TAXONOMY_QZA="${OUTPUT_DIR}/taxonomy-homd.qza"
TAXONOMY_QZV="${OUTPUT_DIR}/taxonomy-homd.qzv"

PHYLOGENY_DIR="${OUTPUT_DIR}/phylogenetic_tree"

mkdir -p "$OUTPUT_DIR"

# ---------------- Sanity checks ----------------

if ! command -v qiime >/dev/null 2>&1; then
  echo "[ERROR] qiime command not found."
  echo "Activate your QIIME 2 conda environment first."
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "[ERROR] Input FASTQ directory does not exist:"
  echo "$INPUT_DIR"
  exit 1
fi

if [[ ! -f "$CLASSIFIER" ]]; then
  echo "[ERROR] HOMD classifier file does not exist:"
  echo "$CLASSIFIER"
  exit 1
fi

echo "========================================================="
echo "QIIME 2 paired-end 16S V3-V4 pipeline"
echo "========================================================="
echo "Input directory:      $INPUT_DIR"
echo "Output directory:     $OUTPUT_DIR"
echo "Classifier:           $CLASSIFIER"
echo "Threads:              $THREADS"
echo "Forward primer:       $FORWARD_PRIMER"
echo "Reverse primer:       $REVERSE_PRIMER"
echo "DADA2 trunc-len-f:    $TRUNC_LEN_F"
echo "DADA2 trunc-len-r:    $TRUNC_LEN_R"
echo "========================================================="

# ==============================================================================
# Step 1: Import FASTQ files to QIIME 2 format
# ==============================================================================

echo "[STEP 1] Importing paired-end FASTQ files..."

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$INPUT_DIR" \
  --input-format 'CasavaOneEightLanelessPerSampleDirFmt' \
  --output-path "$DEMUX_QZA"

echo "[DONE] Imported reads:"
echo "$DEMUX_QZA"

# ==============================================================================
# Step 2: Primer trimming using cutadapt
# ==============================================================================

echo "[STEP 2] Trimming primers with cutadapt..."

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "$DEMUX_QZA" \
  --p-front-f "$FORWARD_PRIMER" \
  --p-front-r "$REVERSE_PRIMER" \
  --p-match-adapter-wildcards \
  --o-trimmed-sequences "$TRIMMED_QZA" \
  --verbose

echo "[DONE] Primer-trimmed reads:"
echo "$TRIMMED_QZA"

# ==============================================================================
# Step 3: Read summary
# ==============================================================================

echo "[STEP 3] Creating demultiplexed read summary..."

qiime demux summarize \
  --i-data "$TRIMMED_QZA" \
  --o-visualization "$TRIMMED_QZV"

echo "[DONE] Read summary:"
echo "$TRIMMED_QZV"

# ==============================================================================
# Step 4: Quality filtering and denoising with DADA2
# ==============================================================================

echo "[STEP 4] Running DADA2 denoise-paired..."

qiime dada2 denoise-paired \
  --p-n-threads "$THREADS" \
  --i-demultiplexed-seqs "$TRIMMED_QZA" \
  --p-trim-left-f "$TRIM_LEFT_F" \
  --p-trim-left-r "$TRIM_LEFT_R" \
  --p-trunc-len-f "$TRUNC_LEN_F" \
  --p-trunc-len-r "$TRUNC_LEN_R" \
  --output-dir "$DADA2_DIR"

echo "[DONE] DADA2 output directory:"
echo "$DADA2_DIR"

# ==============================================================================
# Step 5: DADA2 denoising statistics
# ==============================================================================

echo "[STEP 5] Creating DADA2 denoising statistics summary..."

qiime metadata tabulate \
  --m-input-file "$DADA2_STATS" \
  --o-visualization "$DADA2_STATS_QZV"

echo "[DONE] DADA2 statistics visualization:"
echo "$DADA2_STATS_QZV"

# ==============================================================================
# Step 6: Feature table and representative sequence summaries
# ==============================================================================

echo "[STEP 6] Creating feature table summary..."

qiime feature-table summarize \
  --i-table "$DADA2_TABLE" \
  --o-visualization "$FEATURE_TABLE_QZV"

echo "[DONE] Feature table summary:"
echo "$FEATURE_TABLE_QZV"

echo "[STEP 6] Creating representative sequence summary..."

qiime feature-table tabulate-seqs \
  --i-data "$DADA2_REP_SEQS" \
  --o-visualization "$REP_SEQS_QZV"

echo "[DONE] Representative sequence summary:"
echo "$REP_SEQS_QZV"

# ==============================================================================
# Step 7: Taxonomic classification using HOMD
# ==============================================================================

echo "[STEP 7] Classifying representative sequences using HOMD..."

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$DADA2_REP_SEQS" \
  --o-classification "$TAXONOMY_QZA" \
  --p-n-jobs "$THREADS" \
  --verbose

echo "[DONE] Taxonomy classification:"
echo "$TAXONOMY_QZA"

echo "[STEP 7] Creating taxonomy visualization..."

qiime metadata tabulate \
  --m-input-file "$TAXONOMY_QZA" \
  --o-visualization "$TAXONOMY_QZV"

echo "[DONE] Taxonomy visualization:"
echo "$TAXONOMY_QZV"

# ==============================================================================
# Step 8: Phylogenetic tree
# ==============================================================================

echo "[STEP 8] Building phylogenetic tree..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$DADA2_REP_SEQS" \
  --output-dir "$PHYLOGENY_DIR"

echo "[DONE] Phylogenetic tree output directory:"
echo "$PHYLOGENY_DIR"

# ==============================================================================
# Final summary
# ==============================================================================

echo "========================================================="
echo "QIIME 2 pipeline completed successfully."
echo "========================================================="
echo "Main outputs:"
echo "  Demux reads:              $DEMUX_QZA"
echo "  Primer-trimmed reads:     $TRIMMED_QZA"
echo "  Trimmed read summary:     $TRIMMED_QZV"
echo "  DADA2 table:              $DADA2_TABLE"
echo "  DADA2 representative seqs: $DADA2_REP_SEQS"
echo "  DADA2 stats:              $DADA2_STATS_QZV"
echo "  Feature table summary:    $FEATURE_TABLE_QZV"
echo "  Rep-seqs summary:         $REP_SEQS_QZV"
echo "  Taxonomy:                 $TAXONOMY_QZA"
echo "  Taxonomy summary:         $TAXONOMY_QZV"
echo "  Phylogenetic tree:        $PHYLOGENY_DIR"
echo "========================================================="