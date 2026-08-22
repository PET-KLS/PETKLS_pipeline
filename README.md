# PET-KLS Pipeline

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-2b6cb0.svg)](https://creativecommons.org/licenses/by/4.0/)
![Modality: amyloid PET](https://img.shields.io/badge/modality-amyloid--PET-00a6a6.svg)
![Method: KLS](https://img.shields.io/badge/method-KLS-7c3aed.svg)

Subject-specific amyloid-PET similarity connectomes derived from the complete distribution of voxel-wise standardized uptake value ratio (SUVR) values within cortical parcels.

PET-KLS represents every cortical region by an estimated SUVR probability density and assigns each pair of regions a symmetric distributional-similarity weight. In the accompanying study, the Schaefer-100 atlas yields one symmetric `100 x 100` weighted matrix per participant.

> **Scope.** This repository starts from an already preprocessed SUVR PET image and a parcellation atlas in the same voxel grid. Motion correction, MRI registration, reference-region normalization, scanner-resolution harmonization, partial-volume correction, graph-property calculation, statistical modelling, and visualization are not performed by the KLS computation shown here.

## Processing overview

![PET-KLS processing flowchart](docs/petkls-processing-flowchart.svg)

The editable SVG and a publication-ready PNG are available in [`docs/`](docs/).

## Inputs

| Input | Requirement |
|---|---|
| SUVR PET NIfTI | Subject-level, preprocessed SUVR image. The accompanying study used amyloid-PET SUVR images. |
| Atlas NIfTI | Integer-labelled atlas in the same shape and voxel grid as the PET image. |
| Number of regions | `100` for the Schaefer-100 analysis; labels are expected to be the integers `1,...,R`. |
| Labels file | Required by the current command-line interface for compatibility, but not read by the KLS calculation itself. |

Input images must contain finite numeric values. The supplied implementation retains voxel values satisfying `value > 0` and a region is used only when more than one such voxel remains; it does not separately remove positive infinities or `NaN` values.

## Software requirements

The computation uses Python 3 with `numpy`, `pandas`, `nibabel`, `scipy`, and `tqdm`; `argparse` and `os` are from the Python standard library. Record or pin the exact Python and package versions used for a production run, because this repository code does not encode an environment lock file.

## Exact implemented mathematics

Let `y_(r,k)` denote positive SUVR values in region `r`, with `k = 1,...,n_r` and `n_r > 1`.

### 1. Regional kernel-density estimates

For each region, SciPy's one-dimensional Gaussian KDE estimates

```math
\hat f_r(x)=\frac{1}{n_r h_r}\sum_{k=1}^{n_r}
K\!\left(\frac{x-y_{r,k}}{h_r}\right),
```

where `K` is Gaussian. The implementation calls `scipy.stats.gaussian_kde(..., bw_method=0.2)`. Therefore, `0.2` is a fixed SciPy covariance factor, not Silverman's rule and not an absolute bandwidth in SUVR units. In one dimension, the kernel variance is the regional sample variance multiplied by `0.2^2`.

The KDE is evaluated on the common grid

```math
x_m = \frac{4m}{999},\qquad m=0,\ldots,999,
```

corresponding to 1,000 equally spaced SUVR points from 0 to 4.

### 2. Numerical PDF normalization

The sampled density is normalized over the finite grid using trapezoidal integration:

```math
A_r=\operatorname{trapz}\!\left(\hat f_r(x_m),x_m\right),\qquad
f_r(x_m)=\frac{\hat f_r(x_m)}{A_r}.
```

If `A_r <= 10^-5`, the implementation returns an all-zero density for that region.

Before divergence calculation, every sampled value is clipped below at `epsilon = 10^-10` and the arrays are normalized by their discrete sums:

```math
p_{r,m}=\frac{\max(f_r(x_m),\varepsilon)}
{\sum_{\ell=0}^{999}\max(f_r(x_\ell),\varepsilon)}.
```

### 3. Directed base-2 KL divergence

For regions `i` and `j`, the implemented discrete divergence is

```math
D_2(p_i\|p_j)=\sum_{m=0}^{999}p_{i,m}
\log_2\!\left(\frac{p_{i,m}}{p_{j,m}}\right).
```

### 4. Symmetric divergence

The two directed divergences are added:

```math
J_2(p_i,p_j)=D_2(p_i\|p_j)+D_2(p_j\|p_i).
```

This is the base-2, numerically stabilized discrete form of the Jeffreys divergence. The implementation uses the sum, not the one-half average.

### 5. Kullback-Leibler similarity

The implemented similarity transform is

```math
\operatorname{KLS}_{ij}=\frac{1}{1+J_2(p_i,p_j)}.
```

Consequently, `KLS_ij` is symmetric and lies in `(0,1]` for a valid finite pair: identical sampled distributions give `1`, while increasingly dissimilar distributions approach `0`. A non-finite divergence or a divergence above `10^4` is stored as similarity `0`.

### 6. Matrix assembly

For `R` regions, all unique pairs `i < j` are evaluated and mirrored:

```math
W_{ij}=W_{ji}=\operatorname{KLS}_{ij}.
```

The current output convention leaves the diagonal at zero rather than storing self-similarity. Missing, invalid, or skipped pairs also remain zero; zero must therefore be interpreted as the implementation's unavailable/degenerate-pair sentinel as well as the limiting low-similarity value.

A fuller derivation and line-by-line parameter correspondence are provided in [`docs/mathematical-specification.md`](docs/mathematical-specification.md).

### Known failure modes

The script does not wrap KDE construction or NIfTI indexing in a subject-level error handler. Shape mismatches, non-finite regional samples, or singular/constant regional samples that `gaussian_kde` cannot factorize can therefore stop execution and should be reported as computation failures rather than silently counted as valid connectomes.

## Command-line interface

The computation script accepts the following arguments:

| Argument | Meaning |
|---|---|
| `--pet` | Path to the subject's SUVR PET NIfTI image. |
| `--atlas` | Path to the integer-labelled atlas NIfTI image. |
| `--labels` | Atlas labels file; currently required but not consumed by the calculation. |
| `--n_rois` | Number of atlas regions; default `100`. |
| `--output_dir` | Destination directory. |
| `--tag` | Text inserted into output filenames. |

Generic invocation:

```bash
python path/to/pipeline_script.py \
  --pet sub-001_pet.nii.gz \
  --atlas Schaefer2018_100Parcels_17Networks.nii.gz \
  --labels Schaefer2018_100Parcels_17Networks.csv \
  --n_rois 100 \
  --output_dir derivatives/pet-kls/sub-001 \
  --tag amyloid
```

Use the repository's checked-in pipeline entry point in place of `path/to/pipeline_script.py`.

## Outputs

For each subject, the pipeline writes:

| File suffix | Contents |
|---|---|
| `_connectome.csv` | Symmetric `R x R` KLS matrix. |
| `_pairwise.csv` | Long table with `ROI_1`, `ROI_2`, and `KLS` for successfully evaluated unique pairs. |
| `_pdfs.csv` | Sampled normalized regional PDFs; columns correspond to the 1,000-point SUVR grid. |

## Reproducibility parameters

| Parameter | Implemented value |
|---|---:|
| KDE implementation | `scipy.stats.gaussian_kde` |
| KDE kernel | Gaussian |
| `bw_method` | `0.2` fixed covariance factor |
| SUVR evaluation range | `[0, 4]` |
| Grid size | `1,000` |
| Positive-voxel filter | `SUVR > 0` |
| Minimum regional observations | `n_r > 1` |
| PDF integration | Trapezoidal rule |
| Minimum PDF area | `10^-5` |
| KL clipping constant | `10^-10` |
| Logarithm | Base 2 |
| Symmetrization | Sum of both directed divergences |
| Similarity transform | `1 / (1 + J_2)` |
| Divergence rejection threshold | `J_2 > 10^4` |
| Matrix diagonal | `0` |

## Interpretation

KLS edges quantify similarity between regional amyloid-PET SUVR distributions. They are statistical similarity weights, not anatomical tracts, axonal connections, effective connectivity, or evidence of directed communication. Any graph-theoretical analysis should therefore be described as topology of a PET-signal similarity network.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Until the accompanying article receives its final bibliographic record, cite the repository and the manuscript:

> Monteiro S, et al. *Amyloid PET Similarity Networks in Preclinical Alzheimer's Disease Reveal Early Dynamic Topological Reorganization.* Manuscript under review, 2026.

The approach was modified and adapted from individualized PET distribution-similarity work by Wang and colleagues:

> Wang M, Jiang J, Yan Z, et al. Individual brain metabolic connectome indicator based on Kullback-Leibler Divergence Similarity Estimation predicts progression from mild cognitive impairment to Alzheimer's dementia. *European Journal of Nuclear Medicine and Molecular Imaging.* 2020;47:2753-2764. [https://doi.org/10.1007/s00259-020-04814-x](https://doi.org/10.1007/s00259-020-04814-x)

## License

This repository is licensed under the [Creative Commons Attribution 4.0 International License](LICENSE). Reuse and adaptation are permitted with attribution, a link to the license, and an indication of changes. See the [license statement](LICENSE-STATEMENT.md) for scope, exclusions, and suggested attribution.
