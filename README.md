# 2DTM Postprocessing

A modular Python package for postprocessing 2D template matching results from cryo-EM workflows (e.g., cisTEM), including 2DTM p-value calculation, particle extraction and filtering.

---

## Installation

```bash
git clone https://github.com/kekexinz/2DTM_postprocess_tool.git
cd 2DTM_postprocess_tool
pip install -e .
```

---

## Usage

### `extract-particles`

Extract initial particle peaks from 2DTM search.

```bash
extract-particles \
  --db_file <cistem.db> \
  --tm_job_id 1 \
  --ctf_job_id 1 \
  --pixel_size 1.0 \
  --output <extracted_peaks.star>
```

**Optional parameters:**

| Flag | Description | Default |
|---|---|---|
| `--metric` | `"zscore"` or `"pval"` | `"pval"` |
| `--metric_cutoff` | Threshold for the chosen metric | `8.0` |
| `--threads` | Number of parallel threads | `4` |
| `--local_max_filter` | Image used for `peak_local_max`: `"snr"` or `"zscore"` | `"zscore"` |
| `--min_peak_radius` | `min_distance` in `peak_local_max` (pixels) | `10` |
| `--exclude_borders` | Border exclusion (pixels) to avoid partial particles | `92` |
| `--quadrants` | `1` (first-quadrant only) or `3` (quadrants 1,2,4; recommended for small particles) | `1` |

---

### `filter-particles`

Filter particles based on image thickness and/or angular invariance.

```bash
filter-particles \
  --star_file <extracted_peaks.star> \
  --db_file <cistem.db> \
  --tm_job_id 1 \
  --ctf_job_id 1 \
  --pixel_size 1.0 \
  --output filtered_peaks.star
```

**Optional parameters:**

| Flag | Description |
|---|---|
| `--avg_cutoff_lb` | Lower bound on angular search CC per-pixel average |
| `--sd_cutoff_ub` | Upper bound on angular search CC per-pixel standard deviation |
| `--snr_cutoff_ub` | Upper bound on 2DTM SNR |
| `--filter_by_image_thickness` | Enable filtering by ctffind5 thickness estimate |
| `--thickness_cutoff_lb` | Lower bound on sample thickness (angstroms) |
| `--thickness_cutoff_ub` | Upper bound on sample thickness (angstroms) |
| `--ctf_fitting_score_lb` | Lower bound on CTF fitting score |
| `--ctf_fitting_score_ub` | Upper bound on CTF fitting score |

---

### `measure-template-bias`

Measure the degree of template bias in a 2DTM reconstruction by comparing a full-template reconstruction to an omit-template reconstruction within the omitted region. Python implementation of the bias metric introduced in [Lucas et al. (2023)](https://elifesciences.org/reviewed-preprints/90486v1), adapted for [Zhang et al. (2026)](https://elifesciences.org/reviewed-preprints/109790).

```bash
measure-template-bias \
  --full_recon <recon_full.mrc> \
  --omit_recon <recon_omit.mrc> \
  --templates <full_template.mrc> <omit_template.mrc>
```

Alternatively, provide a precomputed difference map instead of both templates:
```bash
measure-template-bias \
  --full_recon <recon_full.mrc> \
  --omit_recon <recon_omit.mrc> \
  --diff_map <diff_template.mrc>
```

**Optional parameters:**

| Flag | Description | Default |
|---|---|---|
| `--threshold_divisor` | Mask tightness: voxels where `diff > avg_top100 / divisor` are included. Lower = tighter mask (recommended: 3-5) | `10` |
| `--top_n` | Number of top voxels used to compute the threshold | `100` |
| `--save_mask` | Save the binary mask as an MRC file | — |
| `--plot_mask` | Plot central Z/Y/X sections of the mask for visual debugging | `False` |

**Output:**
- Degree of bias: fraction of the full reconstruction's density in the omitted region that is template-dependent (0 = no bias, 1 = fully biased)
- Correlation coefficient between full and omit reconstructions
- Optional: saved binary mask (`.mrc`) and central-section plot (`.png`)

---

### 3D reconstruction & refinement in cisTEM

The output `extracted_peaks.star` and `filtered_peaks.star` can be imported into cisTEM as a RefinementPackage for further 3D reconstruction and refinement.
