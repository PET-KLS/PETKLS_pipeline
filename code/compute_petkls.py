#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.stats import gaussian_kde
from scipy.integrate import trapezoid
import argparse
from tqdm import tqdm

# ------------------------ KLS Utilities ------------------------

def safe_kl_divergence(p, q, eps=1e-10):
    p = np.clip(p, eps, None)
    q = np.clip(q, eps, None)
    p /= np.sum(p)
    q /= np.sum(q)
    return np.sum(p * np.log2(p / q))

def symmetric_kl_divergence(p, q, eps=1e-10):
    try:
        return safe_kl_divergence(p, q, eps) + safe_kl_divergence(q, p, eps)
    except Exception as e:
        print(f"[WARN] Symmetric KL error: {e}")
        return np.nan

def kl_to_similarity(skl):
    if not np.isfinite(skl) or skl > 1e4:
        return 0.0
    return 1 / (1 + skl)

def normalize_pdf(pdf_vals, x_vals):
    area = trapezoid(pdf_vals, x_vals)
    return pdf_vals / area if area > 1e-5 else np.zeros_like(pdf_vals)

# ------------------------ Core Processing ------------------------

def extract_roi_voxels(pet_data, atlas_data, n_rois):
    roi_voxels = {}
    for roi_id in range(1, n_rois + 1):
        mask = (atlas_data == roi_id)
        values = pet_data[mask]
        values = values[values > 0]  # remove background
        if values.size > 1:
            roi_voxels[roi_id] = values
    return roi_voxels

def compute_kls_matrix_and_pdfs(roi_voxels, n_rois, x_vals):
    kls_matrix = np.zeros((n_rois, n_rois))
    records = []
    kde_dict = {}
    pdf_dict = {}

    for roi_id, values in roi_voxels.items():
        kde = gaussian_kde(values, bw_method=0.2)
        pdf_vals = kde(x_vals)
        norm_pdf = normalize_pdf(pdf_vals, x_vals)
        kde_dict[roi_id] = kde
        pdf_dict[roi_id] = norm_pdf

    for i in tqdm(range(1, n_rois + 1), desc="Computing KLS"):
        for j in range(i + 1, n_rois + 1):
            if i in pdf_dict and j in pdf_dict:
                p = pdf_dict[i]
                q = pdf_dict[j]

                if np.all(p < 1e-6) or np.all(q < 1e-6):
                    continue

                skl = symmetric_kl_divergence(p, q)
                kls = kl_to_similarity(skl)

                kls_matrix[i - 1, j - 1] = kls
                kls_matrix[j - 1, i - 1] = kls

                records.append({
                    "ROI_1": i,
                    "ROI_2": j,
                    "KLS": kls
                })

    pairwise_df = pd.DataFrame(records)
    pdf_df = pd.DataFrame.from_dict(pdf_dict, orient="index")
    pdf_df.columns = [f"SUVR_{i}" for i in range(len(x_vals))]
    pdf_df.index.name = "ROI"

    return kls_matrix, pairwise_df, pdf_df

# ------------------------ Main Script ------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute SUVR-based KLS connectome")
    parser.add_argument('--pet', required=True, help='Path to SUVR PET NIfTI file')
    parser.add_argument('--atlas', required=True, help='Path to parcellation atlas NIfTI file')
    parser.add_argument('--labels', required=True, help='Path to atlas label names (not used but required)')
    parser.add_argument('--n_rois', type=int, default=100, help='Number of ROIs')
    parser.add_argument('--output_dir', required=True, help='Directory to store output files')
    parser.add_argument('--tag', required=True, help='Custom tag to include in output filenames')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    x_vals = np.linspace(0, 4, 1000)

    subject_prefix = os.path.basename(args.pet).replace("_pet.nii.gz", "")
    print(f"[INFO] Processing subject: {subject_prefix}")

    pet_data = nib.load(args.pet).get_fdata()
    atlas_data = nib.load(args.atlas).get_fdata()

    roi_voxels = extract_roi_voxels(pet_data, atlas_data, args.n_rois)
    kls_matrix, pairwise_df, pdf_df = compute_kls_matrix_and_pdfs(roi_voxels, args.n_rois, x_vals)

    base = f"{subject_prefix}_{args.tag}_atlas-KLS"

    pd.DataFrame(kls_matrix).to_csv(os.path.join(args.output_dir, f"{base}_connectome.csv"), index=False)
    pairwise_df.to_csv(os.path.join(args.output_dir, f"{base}_pairwise.csv"), index=False)
    pdf_df.to_csv(os.path.join(args.output_dir, f"{base}_pdfs.csv"))

    print(f"[SAVED] Matrix, pairwise, and PDF CSVs for {subject_prefix}")

if __name__ == '__main__':
    main()
