# 2DTM Postprocessing

A modular Python package for postprocessing 2D template matching results from cryo-EM workflows (e.g., cisTEM), including 2DTM p-value calculation, particle extraction and filtering.

---

## Installation

```bash
git clone https://github.com/kekexinz/2DTM_postprocess_tool.git
cd 2DTM_postprocess_tool
pip install -e . # editable mode
```

## 📦 Usage

### `extract-particles`
Extract initial particle peaks from 2DTM search.
```bash
extract-particles \ 
--db_file <cistem.db> \
--tm_job_id 1 \
--ctf_job_id 1 \
--pixel_size 1.0 \
--output <extracted_peaks.star>
[--metric pval] \ # "zscore" or "pval"
[--metric_cutoff 8.0] \
[--threads 22] \
[--local_max_filter] \ # "snr" or "zscore" (default) used for skimage peak_local_max
[--min_peak_radius 10] \ # used for "min_distance" in skimage peak_local_max
[--exclude_borders 92] \ # avoid finding partial particles near the edge of the image, used for skimage peak_local_max
[--quadrants 1] \ # 1 (default) or 3, calculating p-value for only the first-quadrant or quadrant 1,2,4 (recommended for small particles) 

```

### `filter-particles`

Filter particles based on image thickness and/or angular invariance.

```bash
filter-particles \
  --star_file <extracted_peaks.star> \ # output from extract-particles
  --db_file <cistem.db> \
  --tm_job_id 1 \
  --ctf_job_id 1 \
  --pixel_size 1.0 \
  --output filtered_peaks.star \
  [--avg_cutoff_lb] \ # angular search CC per-pixel avg
  [--sd_cutoff_ub] \ # angular search CC per-pixel sd
  [--snr_cutoff_ub] \ 
  [--filter_by_image_thickness] \ # ctffind5 parameters
  [--thickness_cutoff_lb] \
  [--thickness_cutoff_ub] \
  [--ctf_fitting_score_lb] \
  [--ctf_fitting_score_ub] \
```

### `measure-template-bias`

Measure the degree of template bias in a 2DTM reconstruction by comparing a full-template reconstruction to an omit-template reconstruction within the omitted region.

```bash
measure-template-bias \
  --full_recon <recon_full.mrc> \
  --omit_recon <recon_omit.mrc> \
  --templates <full_template.mrc> <omit_template.mrc> \
  [--threshold_divisor 5] \
  [--save_mask mask.mrc] \
  [--plot_mask] \
  [--top_n 100]
```

Alternatively, provide a precomputed difference map instead of both templates:
```bash
measure-template-bias \
  --full_recon <recon_full.mrc> \
  --omit_recon <recon_omit.mrc> \
  --diff_map <diff_template.mrc> \
  [--threshold_divisor 5]
```

**Output:**
- Degree of bias: fraction of the full reconstruction's density in the omitted region that is template-dependent (0 = no bias, 1 = fully biased)
- Correlation coefficient between full and omit reconstructions
- Optional: saved binary mask (`.mrc`) and central-section plot (`.png`) for debugging

**Parameters:**
- `--threshold_divisor`: Controls mask tightness. Mask includes voxels where `diff > avg_top100 / divisor`. Lower values = tighter mask on the omitted atoms (default: 10, recommended: 3-5)
- `--plot_mask`: Plot central Z/Y/X sections of the mask for visual inspection
- `--top_n`: Number of top voxels used to compute the threshold (default: 100)

### 3D reconstruction & refinement in cisTEM
The output extracted_peaks.star and filtered_peaks.star can be imported into cisTEM as a RefinementPackage for further 3D reconstruction and refinement.
