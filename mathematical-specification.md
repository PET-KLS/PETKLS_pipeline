# PET-KLS mathematical and numerical specification

This document describes the supplied implementation exactly. It distinguishes the continuous statistical definitions from the finite-grid operations executed by the code.

## Notation

- `R`: number of atlas regions (`R = 100` in the accompanying Schaefer-100 study).
- `y_(r,k)`: the `k`th retained SUVR voxel in region `r`.
- `n_r`: number of retained voxels in region `r`.
- `M = 1000`: number of evaluation-grid points.
- `x_m = 4m/(M-1)`, for `m = 0,...,M-1`.
- `epsilon = 10^-10`: lower clipping constant used in the divergence calculation.

## 1. Regional voxel extraction

For atlas label `r`, the raw regional vector is

```math
Y_r^{\mathrm{raw}}=\{I(v):A(v)=r\},
```

where `I` is the SUVR image and `A` the integer atlas. The implementation retains

```math
Y_r=\{y\in Y_r^{\mathrm{raw}}:y>0\}.
```

A regional density is attempted only when `n_r = |Y_r| > 1`. The code does not independently resample or register the atlas; PET and atlas arrays must already have compatible grids and shapes.

## 2. Gaussian KDE

For a one-dimensional regional sample, the conceptual KDE is

```math
\hat f_r(x)=\frac{1}{n_r h_r}\sum_{k=1}^{n_r}
\frac{1}{\sqrt{2\pi}}
\exp\!\left[-\frac{1}{2}\left(\frac{x-y_{r,k}}{h_r}\right)^2\right].
```

The call `gaussian_kde(Y_r, bw_method=0.2)` sets SciPy's covariance factor to `c = 0.2`. In one dimension, the kernel covariance is

```math
\Sigma_{K,r}=c^2\Sigma_{Y,r}=0.2^2\Sigma_{Y,r}.
```

Thus, the effective kernel standard deviation is proportional to the within-region sample standard deviation. The scalar is not Silverman's rule and is not a constant absolute SUVR bandwidth.

The KDE is sampled only on `[0,4]`. Probability mass outside this interval is not represented in subsequent calculations.

## 3. Finite-interval and discrete normalization

Let

```math
u_{r,m}=\hat f_r(x_m).
```

First, the array is normalized using the trapezoidal area on `[0,4]`:

```math
A_r=\operatorname{trapz}(u_r,x),\qquad
\bar u_{r,m}=u_{r,m}/A_r,
```

provided `A_r > 10^-5`; otherwise `\bar u_r` is the zero vector.

For KL calculation, the code then constructs discrete probability masses

```math
p_{r,m}=\frac{\max(\bar u_{r,m},\varepsilon)}
{\sum_{\ell=0}^{M-1}\max(\bar u_{r,\ell},\varepsilon)}.
```

Because all grid spacings are equal, this second normalization yields the finite discrete approximation used by the summation. It also guarantees strictly positive denominators after clipping.

## 4. Directed KL divergence in bits

The implementation uses a base-2 logarithm:

```math
D_2(p_i\|p_j)=\sum_{m=0}^{M-1}p_{i,m}
\log_2\left(\frac{p_{i,m}}{p_{j,m}}\right).
```

The divergence is directional: in general,

```math
D_2(p_i\|p_j)\ne D_2(p_j\|p_i).
```

## 5. Symmetric divergence

The implemented symmetric dissimilarity is

```math
J_2(p_i,p_j)=D_2(p_i\|p_j)+D_2(p_j\|p_i).
```

This is a base-2 Jeffreys divergence. No factor of `1/2` is applied. For valid normalized arrays, `J_2 >= 0`, it is symmetric, and it is zero exactly when the sampled distributions are equal.

## 6. Similarity transformation

The dissimilarity is mapped monotonically to similarity by

```math
\operatorname{KLS}_{ij}=\frac{1}{1+J_2(p_i,p_j)}.
```

Properties:

1. `0 < KLS_ij <= 1` for a valid finite divergence.
2. `KLS_ij = KLS_ji`.
3. `KLS_ij = 1` when the discretized distributions are identical.
4. `KLS_ij -> 0` as `J_2 -> infinity`.

The implementation returns `0` when `J_2` is non-finite or exceeds `10^4`.

## 7. Pair eligibility and sentinel zeros

A pair is calculated only if both regional PDF arrays exist and neither array satisfies `all(pdf < 10^-6)`. Otherwise, no long-format pair record is emitted and its matrix entry remains zero.

The initialized matrix also retains a zero diagonal. Therefore:

- diagonal zero means “self-edge not stored”;
- off-diagonal zero can mean invalid/unavailable/skipped pair or rejected divergence;
- a valid finite KLS value is strictly positive.

Downstream code should not automatically interpret all zeros as measured biological dissimilarity without checking pair eligibility.

## 8. Output invariants

For a valid `R`-region result:

```math
W\in\mathbb{R}^{R\times R},\qquad W=W^\top,\qquad
\operatorname{diag}(W)=0,
```

with off-diagonal entries in `[0,1]`. For `R = 100`, there are at most

```math
\binom{100}{2}=4950
```

unique undirected pairs.

## 9. Parameter-to-code correspondence

| Mathematical role | Implemented expression |
|---|---|
| Evaluation grid | `np.linspace(0, 4, 1000)` |
| Regional inclusion | `values = values[values > 0]`; `values.size > 1` |
| KDE factor | `gaussian_kde(values, bw_method=0.2)` |
| Continuous-grid normalization | `trapezoid(pdf_vals, x_vals)`; area floor `1e-5` |
| KL stabilization | `np.clip(..., 1e-10, None)` |
| Discrete normalization | division by array sum |
| Directed divergence | `sum(p * log2(p / q))` |
| Symmetric divergence | directed KL in both directions, added |
| Similarity transform | `1 / (1 + skl)` |
| Divergence guard | invalid or `skl > 1e4` maps to `0` |
| Matrix convention | symmetric off-diagonal values; zero diagonal |

## 10. Scope boundary

This specification covers regional distribution estimation and KLS matrix construction only. It does not specify upstream PET preprocessing or downstream graph construction. In particular, similarity-to-distance conversion, matrix thresholding/normalization, disconnected-path treatment, and graph metrics must be documented from the exact graph-property implementation used for a downstream analysis.

## 11. Reconciliation of the legacy schematic and manuscript description with the code

The corrected repository figure follows the executable implementation. These are the material corrections to the earlier schematic:

| Legacy schematic/manuscript element | Implemented PET-KLS computation |
|---|---|
| FDG-PET/metabolic input | Preprocessed amyloid-PET SUVR image in the accompanying study |
| AAL atlas with 90 regions | Schaefer atlas with 100 cortical parcels in the accompanying study; the script accepts `R` through `--n_rois` |
| MNI normalization, Gaussian smoothing, and SUVR calculation inside the pipeline | The script starts from an already preprocessed SUVR image and performs none of those upstream operations |
| Silverman's bandwidth rule | Fixed SciPy KDE covariance factor `bw_method=0.2` |
| Continuous quadrature expression | KDE evaluation on 1,000 points from 0 to 4, followed by a discrete base-2 KL sum |
| Manuscript equation using one half of the two-direction sum | Sum of the two directed divergences, with no factor of one half |
| Exponential conversion `exp(-J)` | Reciprocal conversion `1/(1+J)` |
| Graph metrics, prediction, and MCE biomarker as pipeline stages | Outside this repository computation; its terminal products are the KLS matrix, pair table, and regional PDF table |

Where narrative descriptions or legacy diagrams differ from these operations, the checked-in executable code and this numerical specification define the reproducible implementation.

## 12. Explicit computation failures

The implementation has no subject-level exception handler around array indexing or `gaussian_kde`. Incompatible PET/atlas shapes, non-finite regional values, or singular regional samples can raise an exception and stop processing. Such cases should be logged and reported as computation failures; they are distinct from a successfully written matrix containing sentinel-zero pairs.
