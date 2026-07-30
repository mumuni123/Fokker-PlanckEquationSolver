#!/usr/bin/env python3
"""Plot separate low-energy and high-energy background-electron spectra.

Unlike EPOCH's two independent ``dist_fn`` diagnostics, this solver writes one
conservative global energy histogram.  The low/high products below are two
windows of that same histogram, so their particle numbers remain directly
comparable and no PIC particle-sampling noise is introduced.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import config
from plot_electron_kinetics import (
    EV_J,
    _safe_time_token,
    _save_csv,
    _select_files,
    _smoothed_spectrum_values,
)
from postprocess_common import save_figure


def _range_mask(
    energy_ev: np.ndarray,
    limits: tuple[float | None, float | None],
) -> np.ndarray:
    lower, upper = limits
    mask = energy_ev > 0.0
    if lower is not None:
        mask &= energy_ev >= float(lower)
    if upper is not None:
        mask &= energy_ev <= float(upper)
    if not np.any(mask):
        raise ValueError(
            f"Energy range {limits} does not overlap "
            f"[{np.min(energy_ev):g}, {np.max(energy_ev):g}] eV"
        )
    return mask


def _distribution_ylabel() -> str:
    mode = str(config.KINETIC_DISTRIBUTION_MODE).lower()
    if mode == "count":
        return "particles per energy bin"
    if mode == "density":
        return r"dN/dE [eV$^{-1}$]"
    if mode == "probability":
        return r"normalized probability density dP/dE [eV$^{-1}$]"
    raise ValueError(
        "KINETIC_DISTRIBUTION_MODE must be 'count', 'density', or 'probability'"
    )


def _plot_window(
    output_dir: Path,
    window_name: str,
    limits: tuple[float | None, float | None],
    selected: list[tuple[Path, float]],
) -> None:
    area_m2 = float(config.KINETIC_TRANSVERSE_AREA_M2)
    if not np.isfinite(area_m2) or area_m2 <= 0.0:
        raise ValueError("KINETIC_TRANSVERSE_AREA_M2 must be finite and positive")

    window_dir = output_dir / f"energy_spectrum_{window_name}"
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    summaries: list[dict[str, object]] = []
    plotted_curve = False

    for index, (path, time_fs) in enumerate(selected, start=1):
        data = np.loadtxt(path, comments="#", ndmin=2)
        if data.shape[1] != 7 or not np.all(np.isfinite(data)):
            raise ValueError(f"{path}: invalid seven-column energy spectrum")

        left, geometric, mean, right, counts_areal, density_areal, fraction = data.T
        widths = right - left
        if np.any(widths <= 0.0) or np.any(counts_areal < 0.0):
            raise ValueError(f"{path}: invalid energy-bin width or negative count")

        counts = counts_areal * area_m2
        density = density_areal * area_m2
        total_number = float(np.sum(counts))
        window_mask = _range_mask(mean, limits)
        window_number = float(np.sum(counts[window_mask]))
        window_energy_j = float(np.sum(counts[window_mask] * mean[window_mask] * EV_J))
        window_fraction = (
            window_number / total_number if total_number > 0.0 else np.nan
        )
        probability_density = (
            counts / total_number / widths
            if total_number > 0.0
            else np.zeros_like(counts)
        )

        plot_y, plotted_counts, _ = _smoothed_spectrum_values(counts, widths)
        plot_mask = window_mask & np.isfinite(plot_y)
        if config.KINETIC_LOG_Y:
            plot_mask &= plot_y > 0.0
            plot_mask &= (
                plotted_counts
                >= config.KINETIC_MIN_BIN_FRACTION * total_number
            )
        if np.any(plot_mask):
            ax.plot(
                mean[plot_mask],
                plot_y[plot_mask],
                linewidth=1.7,
                label=f"t={time_fs:.3g} fs",
            )
            plotted_curve = True

        significant = window_mask & (
            counts >= config.KINETIC_MIN_BIN_FRACTION * total_number
        )
        maximum_resolved_energy = (
            float(np.max(right[significant])) if np.any(significant) else np.nan
        )
        summaries.append(
            {
                "file": path.name,
                "time_fs": time_fs,
                "transverse_area_m2": area_m2,
                "window_lower_eV": limits[0],
                "window_upper_eV": limits[1],
                "total_background_electrons": total_number,
                "window_electrons": window_number,
                "window_particle_fraction": window_fraction,
                "window_kinetic_energy_J": window_energy_j,
                "window_mean_energy_eV": (
                    window_energy_j / window_number / EV_J
                    if window_number > 0.0
                    else np.nan
                ),
                "fraction_above_100keV": (
                    float(np.sum(counts[mean >= 1.0e5])) / total_number
                    if total_number > 0.0
                    else np.nan
                ),
                "fraction_above_1MeV": (
                    float(np.sum(counts[mean >= 1.0e6])) / total_number
                    if total_number > 0.0
                    else np.nan
                ),
                "maximum_resolved_energy_eV": maximum_resolved_energy,
            }
        )

        _save_csv(
            window_dir
            / f"electron_energy_{window_name}_t{_safe_time_token(time_fs)}fs.csv",
            {
                "energy_left_eV": left[window_mask],
                "energy_geometric_center_eV": geometric[window_mask],
                "energy_mean_in_bin_eV": mean[window_mask],
                "energy_right_eV": right[window_mask],
                "N_bin": counts[window_mask],
                "N_bin_m-2": counts_areal[window_mask],
                "dN_dE_eV-1": density[window_mask],
                "dN_dE_m-2_eV-1": density_areal[window_mask],
                "fraction_per_bin": fraction[window_mask],
                "probability_density_eV-1": probability_density[window_mask],
            },
        )
        print(
            f"[{window_name} spectrum {index:2d}/{len(selected)}] "
            f"{path.name}: t={time_fs:g} fs"
        )

    lower, upper = limits
    if lower is not None:
        ax.set_xlim(left=max(float(lower), np.finfo(float).tiny))
    if upper is not None:
        ax.set_xlim(right=float(upper))
    ax.set_xscale("log")
    if config.KINETIC_LOG_Y:
        ax.set_yscale("log")
    ax.set_xlabel("background-electron kinetic energy [eV]")
    ax.set_ylabel(_distribution_ylabel())
    title = "Low-energy background-electron spectrum"
    if window_name == "high":
        title = "High-energy background-electron tail"
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    if plotted_curve:
        ax.legend()
    save_figure(
        fig,
        window_dir / f"electron_energy_spectrum_{window_name}.png",
        config.DPI,
    )
    plt.close(fig)
    _save_csv(
        window_dir / f"electron_energy_spectrum_{window_name}_summary.csv",
        {name: [row[name] for row in summaries] for name in summaries[0]},
    )


def main() -> None:
    selected = _select_files(
        "energy_spectrum_bkg_e_*.dat",
        config.KINETIC_SPECTRUM_FILES,
    )
    output_dir = Path(config.KINETIC_RESULTS_DIR)
    _plot_window(
        output_dir,
        "low",
        tuple(config.KINETIC_LOW_ENERGY_RANGE_EV),
        selected,
    )
    _plot_window(
        output_dir,
        "high",
        tuple(config.KINETIC_HIGH_ENERGY_RANGE_EV),
        selected,
    )
    print("Low/high energy-spectrum diagnostics completed.")


if __name__ == "__main__":
    main()
