#!/usr/bin/env python3
"""Plot one or many normalized background-electron density snapshots."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    figure_output_path,
    filename_token,
    normalize_density_columns,
    read_snapshot_table,
    resolve_input_files,
    save_figure,
    validate_selected_data,
    x_range_mask,
)


def _load_snapshot(path: Path) -> tuple[float | None, list[str], np.ndarray, int]:
    time_fs, labels, data = read_snapshot_table(path)
    try:
        column = labels.index(config.BACKGROUND_DENSITY_COLUMN)
    except ValueError as error:
        raise ValueError(
            f"{path}: column {config.BACKGROUND_DENSITY_COLUMN!r} not found; "
            f"available columns: {', '.join(labels[1:])}"
        ) from error
    if column == 0:
        raise ValueError("The x coordinate cannot be used as a density column")
    labels, data = normalize_density_columns(
        labels,
        data,
        [column],
        config.PARAMETERS_FILE,
    )
    validate_selected_data(path, labels, data, [column])
    return time_fs, labels, data, column


def _time_label(time_fs: float | None, path: Path) -> str:
    return f"t = {time_fs:.3g} fs" if time_fs is not None else path.stem


def _plot_single(
    path: Path,
    time_fs: float | None,
    labels: list[str],
    data: np.ndarray,
    column: int,
) -> None:
    mask, x_limits = x_range_mask(
        data[:, 0],
        config.BACKGROUND_DENSITY_X_AXIS_RANGE,
        "BACKGROUND_DENSITY_X_AXIS_RANGE",
    )
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    ax.plot(data[mask, 0], data[mask, column], linewidth=1.8)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[column])
    ax.set_title(f"Background-electron density, {_time_label(time_fs, path)}")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    output = config.RESULTS_DIR / (
        f"{path.stem}_background_density_{filename_token(labels[column])}.png"
    )
    save_figure(fig, output, config.DPI)
    plt.close(fig)


def _plot_overlay(
    files: list[Path],
    snapshots: list[tuple[float | None, list[str], np.ndarray, int]],
) -> None:
    reference_labels = snapshots[0][1]
    reference_x = snapshots[0][2][:, 0]
    reference_column = snapshots[0][3]
    mask, x_limits = x_range_mask(
        reference_x,
        config.BACKGROUND_DENSITY_X_AXIS_RANGE,
        "BACKGROUND_DENSITY_X_AXIS_RANGE",
    )

    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(snapshots)))
    for color, path, snapshot in zip(colors, files, snapshots):
        time_fs, labels, data, column = snapshot
        if labels != reference_labels or data.shape != snapshots[0][2].shape:
            raise ValueError(f"{path} is incompatible with {files[0]}")
        if not np.allclose(data[:, 0], reference_x):
            raise ValueError(f"{path} uses a different x grid from {files[0]}")
        ax.plot(
            reference_x[mask],
            data[mask, column],
            linewidth=1.6,
            color=color,
            label=_time_label(time_fs, path),
        )

    ax.set_xlabel(reference_labels[0])
    ax.set_ylabel(reference_labels[reference_column])
    ax.set_title("Background-electron density evolution")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(snapshots) <= 12:
        ax.legend(title="snapshot", fontsize="small")
    output = figure_output_path(
        config.RESULTS_DIR,
        files,
        "background_density_evolution",
        filename_token(reference_labels[reference_column]),
    )
    save_figure(fig, output, config.DPI)
    plt.close(fig)


def main() -> None:
    files = resolve_input_files(
        config.BACKGROUND_DENSITY_FILES,
        "BACKGROUND_DENSITY_FILES",
    )
    snapshots = [_load_snapshot(path) for path in files]
    print(f"Processing {len(files)} background-density snapshot(s)")
    for path, snapshot in zip(files, snapshots):
        _plot_single(path, *snapshot)
    if config.BACKGROUND_DENSITY_OVERLAY and len(files) > 1:
        _plot_overlay(files, snapshots)


if __name__ == "__main__":
    main()
