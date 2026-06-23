#!/usr/bin/env python3
"""Plot background plasma density evolution from configured density snapshots."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    figure_output_path,
    filename_token,
    normalize_density_columns,
    parse_header_labels,
    save_figure,
    validate_selected_data,
    x_range_mask,
)


def normalize_density_files(files_config) -> list[Path]:
    if isinstance(files_config, (str, Path)):
        files = [Path(files_config)]
    else:
        files = [Path(path) for path in files_config]

    if not files:
        raise ValueError("BACKGROUND_DENSITY_FILES is empty")
    return sorted(files)


def read_density_file(path: Path) -> tuple[float, list[str], np.ndarray]:
    with path.open("r", encoding="utf-8") as handle:
        time_line = handle.readline().strip()
        header_line = handle.readline().strip()

    match = re.search(r"time\[fs\]\s*=\s*([0-9.+\-eE]+)", time_line)
    time_fs = float(match.group(1)) if match else 0.0

    labels = parse_header_labels(header_line)
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(labels):
        raise ValueError(
            f"{path} has {data.shape[1]} data columns, but header has {len(labels)} labels"
        )
    return time_fs, labels, data


def column_index(labels: list[str], column_name: str) -> int:
    lookup = {label: i for i, label in enumerate(labels)}
    if column_name not in lookup:
        choices = ", ".join(labels[1:])
        raise ValueError(f"Column {column_name!r} not found. Available data columns: {choices}")
    if lookup[column_name] == 0:
        raise ValueError("The x coordinate column cannot be plotted as density")
    return lookup[column_name]


def main() -> None:
    files = normalize_density_files(config.BACKGROUND_DENSITY_FILES)

    snapshots = []
    reference_labels: list[str] | None = None
    reference_x: np.ndarray | None = None
    density_col = -1

    for path in files:
        time_fs, labels, data = read_density_file(path)
        current_col = column_index(labels, config.BACKGROUND_DENSITY_COLUMN)

        # Normalize density by reference density from parameters.h (e.g. n_bkg_e → n_bkg_e / n_e0).
        labels, data = normalize_density_columns(
            labels, data, [current_col], config.PARAMETERS_FILE
        )

        if reference_labels is None:
            reference_labels = labels
            reference_x = data[:, 0]
            density_col = current_col
        else:
            if labels != reference_labels:
                raise ValueError(f"{path} uses different columns from {files[0]}")
            if data.shape[0] != reference_x.shape[0] or not np.allclose(data[:, 0], reference_x):
                raise ValueError(f"{path} uses a different x grid from {files[0]}")

        snapshots.append((time_fs, path, data[:, current_col]))

    if reference_labels is None or reference_x is None:
        raise ValueError("No density snapshots were loaded")

    mask, x_limits = x_range_mask(
        reference_x,
        config.BACKGROUND_DENSITY_X_AXIS_RANGE,
        "BACKGROUND_DENSITY_X_AXIS_RANGE",
    )
    plot_x = reference_x[mask]

    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(snapshots)))
    for color, (time_fs, path, values) in zip(colors, snapshots):
        plot_values = values[mask]
        validate_selected_data(
            path,
            [reference_labels[0], reference_labels[density_col]],
            np.column_stack((plot_x, plot_values)),
            [1],
        )
        ax.plot(plot_x, plot_values, linewidth=1.8, color=color, label=f"{time_fs:.2f} fs")

    ax.set_xlabel(reference_labels[0])
    ax.set_ylabel(reference_labels[density_col])
    ax.set_title("Background plasma density evolution")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(snapshots) <= 12:
        ax.legend(title="time", fontsize="small")

    column_name = filename_token(reference_labels[density_col])
    output = figure_output_path(
        config.RESULTS_DIR,
        files,
        "background_density_evolution",
        column_name,
    )
    save_figure(fig, output, config.DPI)


if __name__ == "__main__":
    main()
