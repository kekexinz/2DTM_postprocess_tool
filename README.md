# 2DTM Postprocessing

A modular Python package for postprocessing 2D template matching results from cryo-EM workflows (e.g., cisTEM), including 2DTM p-value calculation, particle extraction and filtering.

---

## Installation

```bash
git clone https://github.com/kekexinz/2dtm_postprocess
cd 2dtm_postporcess
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
--output <output.star>
[--metric pval] \ # "zscore" or "pval"
[--metric_cutoff 8.0] \
[--threads 22] \
[--local_max_filter] \ # "snr" or "zscore" (default) used for skimage peak_local_max
[--min_peak_radius 10] \ # for finding local max
[--exclude_borders 92] \ # avoid finding partial particles near the edge of the image
[--quadrants 3 ] \ # calculating p-value for only the first quadrant or quadrant 1,2,4 (small proteins) 

```

### `filter-particles`

Filter particles based on image thickness and/or angular invariance.

```bash
filter-particles \
  --star_file <input.star> \
  --db_file <cistem.db> \
  --tm_job_id 1 \
  --ctf_job_id 1 \
  --pixel_size 1.0 \
  --output filtered_output.star \
  [--avg_cutoff_lb 0.0] \ # angular search CC per-pixel avg
  [--sd_cutoff_ub 1.1] \ # angular search CC per-pixel sd
  [--snr_cutoff_ub 6.0] \ 
  [--filter_by_image_thickness] \ # ctffind5 parameters
  [--thickness_cutoff_lb 100.0] \
  [--thickness_cutoff_ub 800.0] \
  [--ctf_fitting_score_lb 0.05] \
  [--ctf_fitting_score_ub 0.2] \
```

### 3D reconstruction & refinement in cisTEM
