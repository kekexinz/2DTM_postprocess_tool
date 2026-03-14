#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# measure_template_bias.py
#
# Python implementation of the template bias metric originally introduced
# in Lucas et al. (2023) and implemented in C++ as MeasureTemplateBiasApp
# in cisTEM. This Python port was developed for the analysis in Zhang
# et al. (2026).
#
# References:
#   Lucas BA et al. (2023) "Baited reconstruction with 2DTM reveals the
#   high-resolution landscape of the cryo-EM structure of the ribosome."
#   eLife. https://elifesciences.org/reviewed-preprints/90486v1
#
#   Zhang K et al. (2026) "Improved cryo-EM reconstruction of sub-50 kDa
#   complexes using 2D template matching."
#   eLife. https://elifesciences.org/reviewed-preprints/109790
# ----------------------------------------------------
from typing import Optional

import argparse
import numpy as np
import mrcfile
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------------- #
# ----------------------------  Helper functions  ----------------------------- #
# ----------------------------------------------------------------------------- #

def read_mrc(path: Path, return_voxel_size: bool = False):
    """
    Read an MRC file into a float32 NumPy array.
    Optionally return voxel_size as well.
    """
    with mrcfile.open(str(path), permissive=True) as mrc:
        data = np.asarray(mrc.data, dtype=np.float32)
        if return_voxel_size:
            return data, mrc.voxel_size.copy()
        return data


def average_of_top_n(vol: np.ndarray, n: int = 100) -> float:
    """
    Return the mean of the largest *n* voxels.
    """
    flat = np.ravel(vol)
    if n >= flat.size:
        return flat.mean()
    return np.mean(np.partition(flat, -n)[-n:])


def correlation(vol1: np.ndarray, vol2: np.ndarray) -> float:
    """
    Pearson correlation coefficient between two volumes.
    """
    v1 = vol1.astype(np.float64).ravel()
    v2 = vol2.astype(np.float64).ravel()
    v1 -= v1.mean()
    v2 -= v2.mean()
    num = np.dot(v1, v2)
    den = np.linalg.norm(v1) * np.linalg.norm(v2)
    return num / den if den != 0 else 0.0


# ----------------------------------------------------------------------------- #
# ---------------------------  Core bias routine  ----------------------------- #
# ----------------------------------------------------------------------------- #

def measure_template_bias(
    full_recon_path: Path,
    omit_recon_path: Path,
    diff_map_path: Optional[Path] = None,
    full_template_path: Optional[Path] = None,
    omit_template_path: Optional[Path] = None,
    top_n: int = 100,
    threshold_divisor: float = 10.0,
    save_mask_path: Optional[Path] = None,
    plot_mask: bool = False,
):
    """
    Compute bias metrics exactly as in the C++ code:

    * mean densities (full / omit / difference) inside a mask of high-value voxels
    * correlation of full & omit reconstructions
    * (optionally) correlation of full & omit templates
    * ratio of those correlations
    * degree of bias
    """

    # ------------------------------------------------------------------ #
    # 1. Load reconstructions
    # ------------------------------------------------------------------ #
    full_recon, voxel_size = read_mrc(full_recon_path, return_voxel_size=True)
    omit_recon = read_mrc(omit_recon_path)

    if full_recon.shape != omit_recon.shape:
        sys.exit("ERROR: Reconstruction volumes do not share dimensions.")

    # ------------------------------------------------------------------ #
    # 2. Decide how to obtain the “difference template” volume
    # ------------------------------------------------------------------ #
    if diff_map_path is not None:
        diff_map = read_mrc(diff_map_path)
        if diff_map.shape != full_recon.shape:
            sys.exit("ERROR: Diff map and recon volumes differ in size.")
        template_corr = None  # we will not compute template correlation
    else:
        if full_template_path is None or omit_template_path is None:
            sys.exit("ERROR: Provide either --diff_map OR both --full_template and --omit_template.")
        full_template = read_mrc(full_template_path)
        omit_template = read_mrc(omit_template_path)
        if full_template.shape != omit_template.shape:
            sys.exit("ERROR: Template volumes differ in size.")
        if full_template.shape != full_recon.shape:
            sys.exit("ERROR: Template and recon volumes differ in size.")
        diff_map = full_template - omit_template
        output_diff_map_path = full_template_path.parent / "diff_from_templates.mrc"
        with mrcfile.new(str(output_diff_map_path), overwrite=True) as mrc:
            mrc.set_data(diff_map.astype(np.float32))
            mrc.voxel_size = voxel_size

        print(f"Saved computed diff map to: {output_diff_map_path}")
        template_corr = correlation(full_template, omit_template)

    # ------------------------------------------------------------------ #
    # 3. Build a mask of high-value voxels in the diff map
    #    (value > 1/10th of the average of top-N voxels)
    # ------------------------------------------------------------------ #
    high_avg = average_of_top_n(diff_map, top_n)
    mask = diff_map > (high_avg / threshold_divisor)
    masked_voxel_count = int(mask.sum())
    if masked_voxel_count == 0:
        sys.exit("ERROR: High-value mask is empty. Check input volumes.")

    if save_mask_path is not None:
        with mrcfile.new(str(save_mask_path), overwrite=True) as mrc:
            mrc.set_data(mask.astype(np.float32))
            mrc.voxel_size = voxel_size
        print(f"Saved mask to: {save_mask_path}")

    if plot_mask:
        mask_float = mask.astype(np.float32)
        mid = [s // 2 for s in mask_float.shape]
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        axes[0].imshow(mask_float[mid[0], :, :], cmap='gray')
        axes[0].set_title(f'Z={mid[0]}')
        axes[1].imshow(mask_float[:, mid[1], :], cmap='gray')
        axes[1].set_title(f'Y={mid[1]}')
        axes[2].imshow(mask_float[:, :, mid[2]], cmap='gray')
        axes[2].set_title(f'X={mid[2]}')
        for ax in axes:
            ax.set_aspect('equal')
        plt.suptitle(f'Mask (divisor={threshold_divisor}, voxels={masked_voxel_count})')
        plt.tight_layout()
        plot_path = save_mask_path.with_suffix('.png') if save_mask_path else Path('mask_plot.png')
        plt.savefig(str(plot_path), dpi=150)
        print(f"Saved mask plot to: {plot_path}")
        plt.close()

    # ------------------------------------------------------------------ #
    # 4. Compute sums / means inside the mask
    # ------------------------------------------------------------------ #
    sum_full = float(full_recon[mask].sum())
    sum_omit = float(omit_recon[mask].sum())
    sum_diff = float((full_recon - omit_recon)[mask].sum())

    mean_full = sum_full / masked_voxel_count
    mean_omit = sum_omit / masked_voxel_count
    mean_difference = sum_diff / masked_voxel_count

    # ------------------------------------------------------------------ #
    # 5. Correlations
    # ------------------------------------------------------------------ #
    recon_corr = correlation(full_recon, omit_recon)

    # ------------------------------------------------------------------ #
    # 6. Assemble and print results
    # ------------------------------------------------------------------ #
    print("\n=== MeasureTemplateBias (Python) ===\n")
    print(f"Masked voxels considered          : {masked_voxel_count}")
    print(f"Mean density (full recon)         : {mean_full: .6g}")
    print(f"Mean density (omit recon)         : {mean_omit: .6g}")
    print(f"Mean density difference           : {mean_difference: .6g}")
    print(f"Correlation (full vs omit recon)  : {recon_corr: .6g}")

    if template_corr is not None:
        ratio_corr = recon_corr / template_corr if template_corr != 0 else np.nan
        print(f"Correlation (full vs omit template) : {template_corr: .6g}")
        print(f"Ratio of recon / template corr      : {ratio_corr: .6g}")

    bias = (sum_full - sum_omit) / sum_full if sum_full != 0 else np.nan
    print(f"\nDegree of bias                    : {bias: .6g}")
    print("\nNormal termination.\n")


# ----------------------------------------------------------------------------- #
# ------------------------------   CLI wrapper   ------------------------------ #
# ----------------------------------------------------------------------------- #

def cli():
    parser = argparse.ArgumentParser(
        description="Re-implementation of MeasureTemplateBiasApp (Relion) in Python.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    required = parser.add_argument_group("required arguments")
    required.add_argument("--full_recon", required=True, type=Path,
                          help="3D reconstruction from full template targets (.mrc)")
    required.add_argument("--omit_recon", required=True, type=Path,
                          help="3D reconstruction from omit template targets (.mrc)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--diff_map", type=Path,
                       help="Difference template map (full − omit) (.mrc)")
    group.add_argument("--templates", nargs=2, metavar=("FULL_TEMPLATE", "OMIT_TEMPLATE"), type=Path,
                       help="Full and omit template volumes (.mrc .mrc)")

    parser.add_argument("--top_n", type=int, default=100,
                        help="Average of top-N voxels used to derive mask threshold")
    parser.add_argument("--threshold_divisor", type=float, default=10.0,
                        help="Voxel is in mask if value > (avg_top_N / divisor)")
    parser.add_argument("--save_mask", type=Path, default=None,
                        help="Save the binary mask as an MRC file")
    parser.add_argument("--plot_mask", action="store_true",
                        help="Plot central sections of the mask for debugging")

    args = parser.parse_args()

    if args.diff_map:
        measure_template_bias(
            full_recon_path=args.full_recon,
            omit_recon_path=args.omit_recon,
            diff_map_path=args.diff_map,
            top_n=args.top_n,
            threshold_divisor=args.threshold_divisor,
            save_mask_path=args.save_mask,
            plot_mask=args.plot_mask,
        )
    else:
        full_template, omit_template = args.templates
        measure_template_bias(
            full_recon_path=args.full_recon,
            omit_recon_path=args.omit_recon,
            full_template_path=full_template,
            omit_template_path=omit_template,
            top_n=args.top_n,
            threshold_divisor=args.threshold_divisor,
            save_mask_path=args.save_mask,
            plot_mask=args.plot_mask,
        )


if __name__ == "__main__":
    cli()
