#!/usr/bin/env python3
"""Background-electron energy, spectrum, and momentum postprocessing.

This is the distribution-function analogue of the EPOCH diagnostics: it reads
the solver's globally integrated moments and deterministic velocity-space
histograms rather than SDF particle samples.  The solver's one-dimensional
areal quantities are converted to totals with KINETIC_TRANSVERSE_AREA_M2.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import read_table, save_figure


EV_J = 1.60218e-19


def _column(labels: list[str], name: str) -> int:
    try:
        return labels.index(name)
    except ValueError as error:
        raise ValueError(
            f"{name!r} is missing; available columns: {', '.join(labels)}"
        ) from error


def _optional_column(
    labels: list[str],
    data: np.ndarray,
    name: str,
) -> np.ndarray:
    if name not in labels:
        return np.full(data.shape[0], np.nan)
    return data[:, labels.index(name)]


def _validate_finite(path: Path, name: str, values: np.ndarray) -> None:
    if not np.all(np.isfinite(values)):
        count = int(np.count_nonzero(~np.isfinite(values)))
        raise ValueError(f"{path}: {name} contains {count} non-finite values")


def _deduplicate_time_rows(time: np.ndarray, data: np.ndarray) -> np.ndarray:
    """Sort by time and retain the last record when a time was written twice."""
    order = np.argsort(time, kind="stable")
    sorted_data = data[order]
    sorted_time = time[order]
    keep = np.ones(sorted_time.size, dtype=bool)
    keep[:-1] = sorted_time[:-1] != sorted_time[1:]
    return sorted_data[keep]


def _save_csv(path: Path, columns: dict[str, Iterable[object]]) -> None:
    names = list(columns)
    arrays = [list(columns[name]) for name in names]
    row_count = len(arrays[0])
    if any(len(array) != row_count for array in arrays):
        raise ValueError(f"CSV columns have inconsistent lengths for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        writer.writerows(zip(*arrays))
    print(f"Saved {path}")


def _plot_curves(
    time_fs: np.ndarray,
    curves: list[tuple[str, np.ndarray]],
    ylabel: str,
    title: str,
    output: Path,
) -> None:
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    for label, values in curves:
        if np.any(np.isfinite(values)):
            ax.plot(time_fs, values, linewidth=1.8, label=label)
    ax.set_xlabel("time [fs]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
    ax.grid(True, alpha=0.3)
    ax.legend()
    save_figure(fig, output, config.DPI)
    plt.close(fig)


def plot_energy_history(output_dir: Path) -> None:
    path = config.DATA_DIR / "scalars.dat"
    labels, data = read_table(path)
    data = _deduplicate_time_rows(data[:, _column(labels, "time[fs]")], data)

    step = data[:, _column(labels, "step")]
    time_fs = data[:, _column(labels, "time[fs]")]
    area_m2 = float(config.KINETIC_TRANSVERSE_AREA_M2)
    if not np.isfinite(area_m2) or area_m2 <= 0.0:
        raise ValueError("KINETIC_TRANSVERSE_AREA_M2 must be finite and positive")

    number_bkg_areal = data[:, _column(labels, "N_bkg_e")]
    ke_bkg_areal = data[:, _column(labels, "KE_bkg_e[J/m2]")]
    number_beam_areal = data[:, _column(labels, "N_beam")]
    ke_beam_areal = data[:, _column(labels, "KE_beam[J/m2]")]
    field_areal = data[:, _column(labels, "E_field[J/m2]")]
    total_areal = data[:, _column(labels, "E_total[J/m2]")]
    injected_areal = _optional_column(labels, data, "E_beam_injected_cum[J/m2]")
    outflow_areal = _optional_column(labels, data, "E_beam_outflow_cum[J/m2]")
    collision_areal = _optional_column(labels, data, "E_collision_cum[J/m2]")
    accounted_areal = _optional_column(labels, data, "E_accounted[J/m2]")
    balance_error_areal = _optional_column(labels, data, "E_balance_error[J/m2]")

    number_bkg = number_bkg_areal * area_m2
    number_beam = number_beam_areal * area_m2
    ke_bkg = ke_bkg_areal * area_m2
    ke_beam = ke_beam_areal * area_m2
    field = field_areal * area_m2
    total = total_areal * area_m2
    injected = injected_areal * area_m2
    outflow = outflow_areal * area_m2
    collision = collision_areal * area_m2
    accounted = accounted_areal * area_m2
    balance_error = balance_error_areal * area_m2

    for name, values in (
        ("time", time_fs),
        ("N_bkg_e", number_bkg),
        ("KE_bkg_e", ke_bkg),
        ("KE_beam", ke_beam),
        ("E_field", field),
        ("E_total", total),
    ):
        _validate_finite(path, name, values)
    if np.any(number_bkg <= 0.0):
        raise ValueError(f"{path}: N_bkg_e must remain positive")

    particle_ke = ke_bkg + ke_beam
    bkg_plus_field = ke_bkg + field
    mean_ke_j = ke_bkg / number_bkg
    mean_ke_ev = mean_ke_j / EV_J
    delta_total = total - total[0]
    total_scale = max(abs(float(total[0])), np.finfo(float).tiny)
    relative_total_change = delta_total / total_scale
    relative_balance_error = balance_error / np.maximum(np.abs(total), total_scale)

    _save_csv(
        output_dir / "energy_history.csv",
        {
            "step": step.astype(np.int64),
            "time_fs": time_fs,
            "transverse_area_m2": np.full(time_fs.shape, area_m2),
            "N_bkg_e": number_bkg,
            "N_bkg_e_m-2": number_bkg_areal,
            "N_beam": number_beam,
            "N_beam_m-2": number_beam_areal,
            "KE_bkg_e_J": ke_bkg,
            "KE_bkg_e_J_m-2": ke_bkg_areal,
            "mean_KE_bkg_e_J": mean_ke_j,
            "mean_KE_bkg_e_eV": mean_ke_ev,
            "KE_beam_J": ke_beam,
            "particle_KE_J": particle_ke,
            "electric_field_energy_J": field,
            "bkg_plus_field_energy_J": bkg_plus_field,
            "total_energy_including_beam_J": total,
            "delta_total_energy_J": delta_total,
            "relative_total_energy_change": relative_total_change,
            "beam_injected_energy_cumulative_J": injected,
            "beam_outflow_energy_cumulative_J": outflow,
            "collision_energy_cumulative_J": collision,
            "accounted_energy_J": accounted,
            "energy_balance_error_J": balance_error,
            "relative_energy_balance_error": relative_balance_error,
        },
    )

    _plot_curves(
        time_fs,
        [
            ("Total energy (including beam)", total),
            ("Total particle kinetic energy", particle_ke),
            ("Electric-field energy", field),
        ],
        "energy [J]",
        "Domain energy history",
        output_dir / "domain_energy_history.png",
    )
    _plot_curves(
        time_fs,
        [("Electric-field energy", field)],
        "electric-field energy [J]",
        "Electric-field energy history",
        output_dir / "electric_field_energy_history.png",
    )
    _plot_curves(
        time_fs,
        [("Background-electron kinetic energy", ke_bkg)],
        "background-electron kinetic energy [J]",
        "Background-electron kinetic-energy history",
        output_dir / "bkgelectron_total_energy_history.png",
    )
    _plot_curves(
        time_fs,
        [("Background-electron mean kinetic energy", mean_ke_ev)],
        "mean kinetic energy [eV]",
        "Background-electron mean kinetic-energy history",
        output_dir / "bkgelectron_mean_energy_history.png",
    )

    fig, (ax_change, ax_error) = plt.subplots(
        2, 1, sharex=True, figsize=(config.FIGSIZE[0], config.FIGSIZE[1] * 1.35)
    )
    ax_change.plot(time_fs, relative_total_change, label=r"$(E_{\rm total}-E_0)/|E_0|$")
    ax_change.set_ylabel("relative change")
    ax_change.grid(True, alpha=0.3)
    ax_change.legend()
    ax_error.plot(time_fs, balance_error, color="tab:red", label="accounting error")
    ax_error.set_xlabel("time [fs]")
    ax_error.set_ylabel("error [J]")
    ax_error.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
    ax_error.grid(True, alpha=0.3)
    ax_error.legend()
    fig.suptitle("Energy change and solver energy accounting")
    save_figure(fig, output_dir / "energy_balance_history.png", config.DPI)
    plt.close(fig)

    print(
        "Energy summary: "
        f"t={time_fs[0]:g}..{time_fs[-1]:g} fs, "
        f"<KE>_bkg={mean_ke_ev[0]:.6g}..{mean_ke_ev[-1]:.6g} eV, "
        f"max |relative balance error|="
        f"{np.nanmax(np.abs(relative_balance_error)):.6e}"
    )


def _time_from_header(path: Path) -> float:
    with path.open("r", encoding="utf-8") as handle:
        for _ in range(8):
            line = handle.readline()
            if not line:
                break
            match = re.search(r"\btime_fs\s+([-+0-9.eE]+)", line)
            if match:
                return float(match.group(1))
    raise ValueError(f"Cannot read time_fs from {path}")


def _metadata_from_header(path: Path) -> dict[str, float]:
    metadata: dict[str, float] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            match = re.match(r"#\s*([A-Za-z0-9_.-]+)\s+([-+0-9.eE]+)\s*$", line)
            if match:
                metadata[match.group(1)] = float(match.group(2))
    return metadata


def _configured_or_available(
    pattern: str,
    configured: list[Path] | None,
) -> list[Path]:
    files = (
        sorted(config.DATA_DIR.glob(pattern))
        if configured is None
        else [Path(item) for item in configured]
    )
    if not files:
        raise FileNotFoundError(f"No files matching {pattern!r} in {config.DATA_DIR}")
    for path in files:
        if not path.exists():
            raise FileNotFoundError(path)
    return files


def _select_files(
    pattern: str,
    configured: list[Path] | None,
) -> list[tuple[Path, float]]:
    files = _configured_or_available(pattern, configured)
    timed = [(path, _time_from_header(path)) for path in files]
    timed.sort(key=lambda item: item[1])
    if configured is not None:
        return timed

    targets = getattr(config, "KINETIC_SAMPLE_TIMES_FS", None)
    if targets is None:
        return [timed[0], timed[-1]] if len(timed) > 1 else timed

    selected: list[tuple[Path, float]] = []
    used: set[Path] = set()
    for target in targets:
        nearest = min(timed, key=lambda item: abs(item[1] - float(target)))
        if nearest[0] not in used:
            selected.append(nearest)
            used.add(nearest[0])
    selected.sort(key=lambda item: item[1])
    return selected


def _mode_values(
    counts: np.ndarray,
    widths: np.ndarray,
    density: np.ndarray,
    coordinate_name: str,
) -> tuple[np.ndarray, str]:
    mode = str(config.KINETIC_DISTRIBUTION_MODE).lower()
    if mode == "count":
        return counts, "particles per bin"
    if mode == "density":
        return density, rf"dN/d{coordinate_name}"
    if mode == "probability":
        total = float(np.sum(counts))
        values = counts / total / widths if total > 0.0 else np.zeros_like(counts)
        return values, rf"normalized probability density dP/d{coordinate_name}"
    raise ValueError(
        "KINETIC_DISTRIBUTION_MODE must be 'count', 'density', or 'probability'"
    )


def _safe_time_token(time_fs: float) -> str:
    return f"{time_fs:.6f}".replace("-", "m").replace(".", "p")


def _smoothed_spectrum_values(
    counts: np.ndarray,
    widths: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Return plotting values using conservative neighboring-bin averaging."""
    window = int(config.KINETIC_SPECTRUM_SMOOTHING_BINS)
    if window < 1 or window % 2 == 0:
        raise ValueError("KINETIC_SPECTRUM_SMOOTHING_BINS must be a positive odd integer")
    kernel = np.ones(window)
    summed_counts = np.convolve(counts, kernel, mode="same")
    summed_widths = np.convolve(widths, kernel, mode="same")
    local_bin_count = np.convolve(np.ones_like(counts), kernel, mode="same")
    mode = str(config.KINETIC_DISTRIBUTION_MODE).lower()
    if mode == "count":
        return summed_counts / local_bin_count, summed_counts, "particles per bin"
    smoothed_density = summed_counts / summed_widths
    if mode == "density":
        return smoothed_density, summed_counts, r"dN/dE [eV$^{-1}$]"
    if mode == "probability":
        total = float(np.sum(counts))
        values = smoothed_density / total if total > 0.0 else np.zeros_like(counts)
        return values, summed_counts, r"normalized probability density dP/dE [eV$^{-1}$]"
    raise ValueError(
        "KINETIC_DISTRIBUTION_MODE must be 'count', 'density', or 'probability'"
    )


def plot_energy_spectra(output_dir: Path) -> None:
    area_m2 = float(config.KINETIC_TRANSVERSE_AREA_M2)
    if not np.isfinite(area_m2) or area_m2 <= 0.0:
        raise ValueError("KINETIC_TRANSVERSE_AREA_M2 must be finite and positive")
    selected = _select_files(
        "energy_spectrum_bkg_e_*.dat",
        config.KINETIC_SPECTRUM_FILES,
    )
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    summaries: list[dict[str, object]] = []

    for index, (path, time_fs) in enumerate(selected, start=1):
        data = np.loadtxt(path, comments="#", ndmin=2)
        if data.shape[1] != 7:
            raise ValueError(f"{path}: expected 7 spectrum columns, got {data.shape[1]}")
        if not np.all(np.isfinite(data)):
            raise ValueError(f"{path}: spectrum contains non-finite values")

        left, geometric, mean, right, counts_areal, density_areal, fraction = data.T
        widths = right - left
        if np.any(widths <= 0.0) or np.any(counts_areal < 0.0):
            raise ValueError(f"{path}: invalid energy-bin width or negative count")
        counts = counts_areal * area_m2
        density = density_areal * area_m2
        y, plotted_counts, ylabel = _smoothed_spectrum_values(counts, widths)
        mask = np.isfinite(y) & (mean > 0.0)
        if config.KINETIC_LOG_Y:
            mask &= y > 0.0
            mask &= (
                plotted_counts
                >= config.KINETIC_MIN_BIN_FRACTION * np.sum(plotted_counts)
            )
        ax.plot(mean[mask], y[mask], linewidth=1.7, label=f"t={time_fs:.3g} fs")

        total = float(np.sum(counts))
        integrated_ke = float(np.sum(counts * mean * EV_J))
        probability_density = (
            counts / total / widths if total > 0.0 else np.zeros_like(counts)
        )
        metadata = _metadata_from_header(path)
        header_total = metadata.get("integrated_N_m-2", np.nan) * area_m2
        header_ke = (
            metadata.get("integrated_exact_kinetic_energy_J_m-2", np.nan)
            * area_m2
        )
        summaries.append(
            {
                "file": path.name,
                "time_fs": time_fs,
                "transverse_area_m2": area_m2,
                "integrated_N": total,
                "header_integrated_N": header_total,
                "relative_N_closure_error": (
                    (total - header_total) / header_total
                    if np.isfinite(header_total) and header_total != 0.0
                    else np.nan
                ),
                "mean_kinetic_energy_eV": integrated_ke / total / EV_J,
                "binned_kinetic_energy_J": integrated_ke,
                "header_exact_kinetic_energy_J": header_ke,
                "relative_energy_binning_error": (
                    (integrated_ke - header_ke) / header_ke
                    if np.isfinite(header_ke) and header_ke != 0.0
                    else np.nan
                ),
            }
        )
        _save_csv(
            output_dir / f"energy_spectrum_t{_safe_time_token(time_fs)}fs.csv",
            {
                "energy_left_eV": left,
                "energy_geometric_center_eV": geometric,
                "energy_mean_in_bin_eV": mean,
                "energy_right_eV": right,
                "N_bin": counts,
                "N_bin_m-2": counts_areal,
                "dN_dE_eV-1": density,
                "dN_dE_m-2_eV-1": density_areal,
                "fraction_per_bin": fraction,
                "probability_density_eV-1": probability_density,
            },
        )
        print(f"[spectrum {index:2d}/{len(selected)}] {path.name}: t={time_fs:g} fs")

    ax.set_xlabel("background-electron kinetic energy [eV]")
    ax.set_ylabel(ylabel)
    ax.set_title("Background-electron energy spectrum")
    ax.set_xscale("log")
    if config.KINETIC_LOG_Y:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    save_figure(fig, output_dir / "electron_energy_spectrum.png", config.DPI)
    plt.close(fig)
    _save_csv(
        output_dir / "energy_spectrum_summary.csv",
        {name: [row[name] for row in summaries] for name in summaries[0]},
    )


def plot_momentum_distributions(output_dir: Path) -> None:
    area_m2 = float(config.KINETIC_TRANSVERSE_AREA_M2)
    if not np.isfinite(area_m2) or area_m2 <= 0.0:
        raise ValueError("KINETIC_TRANSVERSE_AREA_M2 must be finite and positive")
    components = (
        (
            "parallel",
            "parallel",
            config.KINETIC_PARALLEL_MOMENTUM_FILES,
        ),
        (
            "perpendicular",
            "perpendicular",
            config.KINETIC_PERPENDICULAR_MOMENTUM_FILES,
        ),
    )
    summaries: list[dict[str, object]] = []

    for component, display, configured in components:
        selected = _select_files(f"momentum_{component}_bkg_e_*.dat", configured)
        fig, ax = plt.subplots(figsize=config.FIGSIZE)
        for index, (path, time_fs) in enumerate(selected, start=1):
            data = np.loadtxt(path, comments="#", ndmin=2)
            if data.shape[1] != 7:
                raise ValueError(
                    f"{path}: expected 7 momentum columns, got {data.shape[1]}"
                )
            if not np.all(np.isfinite(data)):
                raise ValueError(f"{path}: momentum distribution has non-finite values")

            (
                left,
                center,
                right,
                p_center,
                counts_areal,
                density_u_areal,
                density_p_areal,
            ) = data.T
            widths = right - left
            if np.any(widths <= 0.0) or np.any(counts_areal < 0.0):
                raise ValueError(f"{path}: invalid momentum-bin width or negative count")
            counts = counts_areal * area_m2
            density_u = density_u_areal * area_m2
            density_p = density_p_areal * area_m2
            y, ylabel = _mode_values(counts, widths, density_u, "u")
            mask = np.isfinite(y)
            if config.KINETIC_LOG_Y:
                mask &= y > 0.0
                mask &= counts >= config.KINETIC_MIN_BIN_FRACTION * np.sum(counts)
            ax.plot(center[mask], y[mask], linewidth=1.7, label=f"t={time_fs:.3g} fs")

            total = float(np.sum(counts))
            mean_u = float(np.sum(counts * center) / total)
            rms_u = float(np.sqrt(np.sum(counts * center * center) / total))
            probability_density = (
                counts / total / widths if total > 0.0 else np.zeros_like(counts)
            )
            summaries.append(
                {
                    "component": component,
                    "file": path.name,
                    "time_fs": time_fs,
                    "transverse_area_m2": area_m2,
                    "integrated_N": total,
                    "mean_u": mean_u,
                    "rms_u": rms_u,
                    "mean_p_kg_m_s-1": float(np.sum(counts * p_center) / total),
                }
            )
            _save_csv(
                output_dir
                / f"momentum_{component}_t{_safe_time_token(time_fs)}fs.csv",
                {
                    "u_left": left,
                    "u_center": center,
                    "u_right": right,
                    "p_center_kg_m_s-1": p_center,
                    "N_bin": counts,
                    "N_bin_m-2": counts_areal,
                    "dN_du": density_u,
                    "dN_du_m-2": density_u_areal,
                    "dN_dp_per_kg_m_s": density_p,
                    "dN_dp_m-2_per_kg_m_s": density_p_areal,
                    "probability_density_du-1": probability_density,
                },
            )
            print(
                f"[momentum {component} {index:2d}/{len(selected)}] "
                f"{path.name}: t={time_fs:g} fs"
            )

        coordinate = r"\parallel" if component == "parallel" else r"\perp"
        ax.set_xlabel(rf"$u_{{{coordinate}}}=p_{{{coordinate}}}/(m_e c)$")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Background-electron {display} momentum distribution")
        if config.KINETIC_LOG_Y:
            ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        save_figure(
            fig,
            output_dir / f"electron_momentum_{component}.png",
            config.DPI,
        )
        plt.close(fig)

    _save_csv(
        output_dir / "momentum_summary.csv",
        {name: [row[name] for row in summaries] for name in summaries[0]},
    )


def main() -> None:
    output_dir = Path(config.KINETIC_RESULTS_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Input:  {config.DATA_DIR.resolve()}")
    print(f"Output: {output_dir.resolve()}")

    plot_energy_history(output_dir)
    tasks = (
        ("energy spectrum", plot_energy_spectra),
        ("momentum distribution", plot_momentum_distributions),
    )
    for name, function in tasks:
        try:
            function(output_dir)
        except FileNotFoundError as error:
            print(f"Skipped {name}: {error}")

    print("\nElectron kinetic diagnostics completed.")
    print("- energy_history.csv and five energy-history figures")
    print("- per-snapshot spectrum/momentum CSV files")
    print("- energy_spectrum_summary.csv and momentum_summary.csv")


if __name__ == "__main__":
    main()
