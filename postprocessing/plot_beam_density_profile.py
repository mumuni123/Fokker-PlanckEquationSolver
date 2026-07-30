#!/usr/bin/env python3
"""Plot one or many normalized beam-electron density snapshots."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    figure_output_path,
    filename_token,
    normalize_columns,
    normalize_density_columns,
    read_snapshot_table,
    resolve_input_files,
    save_figure,
    validate_selected_data,
    x_range_mask,
    y_axis_label,
)


def _load_snapshot(path: Path) -> tuple[float | None, list[str], np.ndarray, list[int]]:
    time_fs, labels, data = read_snapshot_table(path)
    columns = normalize_columns(labels, config.DENSITY_COLUMNS)
    labels, data = normalize_density_columns(
        labels,
        data,
        columns,
        config.PARAMETERS_FILE,
    )
    validate_selected_data(path, labels, data, columns)
    return time_fs, labels, data, columns


def _time_label(time_fs: float | None, path: Path) -> str:
    return f"t = {time_fs:.3g} fs" if time_fs is not None else path.stem


def _plot_single(
    path: Path,
    time_fs: float | None,
    labels: list[str],
    data: np.ndarray,
    columns: list[int],
) -> None:
    mask, x_limits = x_range_mask(
        data[:, 0],
        config.DENSITY_X_AXIS_RANGE,
        "DENSITY_X_AXIS_RANGE",
    )
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    for column in columns:
        ax.plot(
            data[mask, 0],
            data[mask, column],
            linewidth=1.8,
            label=labels[column],
        )
    ax.set_xlabel(labels[0])
    ax.set_ylabel(y_axis_label(labels, columns))
    ax.set_title(f"Beam-electron density, {_time_label(time_fs, path)}")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(columns) > 1:
        ax.legend()
    token = "_".join(filename_token(labels[column]) for column in columns)
    save_figure(
        fig,
        config.RESULTS_DIR / f"{path.stem}_beam_density_{token}.png",
        config.DPI,
    )
    plt.close(fig)


def _plot_overlay(
    files: list[Path],
    snapshots: list[tuple[float | None, list[str], np.ndarray, list[int]]],
) -> None:
    reference_labels = snapshots[0][1]
    reference_data = snapshots[0][2]
    reference_columns = snapshots[0][3]
    if len(reference_columns) != 1:
        raise ValueError(
            "BEAM_DENSITY_OVERLAY requires exactly one DENSITY_COLUMNS entry"
        )
    mask, x_limits = x_range_mask(
        reference_data[:, 0],
        config.DENSITY_X_AXIS_RANGE,
        "DENSITY_X_AXIS_RANGE",
    )
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    colors = plt.cm.plasma(np.linspace(0.0, 1.0, len(snapshots)))
    for color, path, (time_fs, labels, data, columns) in zip(
        colors, files, snapshots
    ):
        if labels != reference_labels or data.shape != reference_data.shape:
            raise ValueError(f"{path} is incompatible with {files[0]}")
        if not np.allclose(data[:, 0], reference_data[:, 0]):
            raise ValueError(f"{path} uses a different x grid from {files[0]}")
        ax.plot(
            data[mask, 0],
            data[mask, columns[0]],
            linewidth=1.6,
            color=color,
            label=_time_label(time_fs, path),
        )
    ax.set_xlabel(reference_labels[0])
    ax.set_ylabel(reference_labels[reference_columns[0]])
    ax.set_title("Beam-electron density evolution")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(snapshots) <= 12:
        ax.legend(title="snapshot", fontsize="small")
    output = figure_output_path(
        config.RESULTS_DIR,
        files,
        "beam_density_evolution",
        filename_token(reference_labels[reference_columns[0]]),
    )
    save_figure(fig, output, config.DPI)
    plt.close(fig)


def main() -> None:
    files = resolve_input_files(config.BEAM_DENSITY_FILES, "BEAM_DENSITY_FILES")
    snapshots = [_load_snapshot(path) for path in files]
    print(f"Processing {len(files)} beam-density snapshot(s)")
    for path, snapshot in zip(files, snapshots):
        _plot_single(path, *snapshot)
    if config.BEAM_DENSITY_OVERLAY and len(files) > 1:
        _plot_overlay(files, snapshots)


if __name__ == "__main__":
    main()
