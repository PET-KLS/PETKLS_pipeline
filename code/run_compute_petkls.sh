#!/usr/bin/env bash
#SBATCH --job-name=indiv_pet_connectome
#SBATCH --output=logs/pet_kls_%j.out
#SBATCH --error=logs/pet_kls_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=16:00:00
#SBATCH --partition=<LONG_PARTITION>
#SBATCH --nice=0

set -euo pipefail

module load Anaconda3
conda activate pet_kls

# ======== CONFIGURATION ======== #
# Replace this placeholder or export PROJECT_ROOT before submission:
# export PROJECT_ROOT="/path/to/project"
PROJECT_ROOT="${PROJECT_ROOT:-/path/to/project}"

DERIVATIVES_DIR="${PROJECT_ROOT}/derivatives"
SUVR_DIR="${DERIVATIVES_DIR}/suvr"
ATLAS_DIR="${PROJECT_ROOT}/atlases/Schaefer2018"
ATLAS="${ATLAS_DIR}/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz"
LABELS="${ATLAS_DIR}/centroids/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv"
OUTPUT_DIR="${DERIVATIVES_DIR}/pet_connectomes_kls"
SCRIPT_PATH="${PROJECT_ROOT}/scripts/compute_indiv_connectome.py"

TAG="desc-suvr_KLS"
N_ROIS=100

mkdir -p logs "$OUTPUT_DIR"

# ======== VALIDATION ======== #
for REQUIRED_PATH in "$SUVR_DIR" "$ATLAS" "$LABELS" "$SCRIPT_PATH"; do
    if [[ ! -e "$REQUIRED_PATH" ]]; then
        echo "[ERROR] Required path does not exist: $REQUIRED_PATH" >&2
        exit 1
    fi
done

# ======== PROCESS PET IMAGES ======== #
find "$SUVR_DIR" -type f -name "*desc-suvrcereb_pet.nii.gz" -print0 |
while IFS= read -r -d '' PET_FILE; do
    BASENAME=$(basename "$PET_FILE")
    SUBJ_ID=$(printf '%s\n' "$BASENAME" | cut -d '_' -f 1)
    SES=$(printf '%s\n' "$BASENAME" | cut -d '_' -f 2)

    SUBJECT_OUTPUT_DIR="${OUTPUT_DIR}/${SUBJ_ID}/${SES}"
    OUT_FILE="${SUBJECT_OUTPUT_DIR}/${SUBJ_ID}_${TAG}_atlas-KLS_connectome.csv"

    mkdir -p "$SUBJECT_OUTPUT_DIR"

    if [[ -f "$OUT_FILE" ]]; then
        echo "[SKIP] Output exists for ${SUBJ_ID}, session ${SES}"
        continue
    fi

    echo "[RUN] Computing connectome for ${SUBJ_ID}, session ${SES}"

    python "$SCRIPT_PATH" \
        --pet "$PET_FILE" \
        --atlas "$ATLAS" \
        --labels "$LABELS" \
        --n_rois "$N_ROIS" \
        --output_dir "$SUBJECT_OUTPUT_DIR" \
        --tag "$TAG"
done

echo "[DONE] PET connectome processing completed."
