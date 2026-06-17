#!/usr/bin/env python3
"""Plot diagnostics from face-centered current snapshots."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    filename_token,
    parse_header_labels,
    save_figure,
    validate_selected_data,
    x_range_mask,
)


def configured_files() -> list[Path]:
    files_config = config.FACE_CURRENT_FILES
    if files_config is None:
        raise ValueError("FACE_CURRENT_FILES must specify one file or a non-empty file list")
    if isinstance(files_config, (str, Path)):
        files = [Path(files_config)]
    else:
        files = [Path(path) for path in files_config]

    if not files:
        raise ValueError("FACE_CURRENT_FILES is empty")
    return sorted(files)


def read_face_current(path: Path) -> tuple[float, list[str], np.ndarray]:
    with path.open("r", encoding="utf-8") as handle:
        time_line = handle.readline().strip()
        header_line = handle.readline().strip()

    match = re.search(r"time\[fs\]\s*=\s*([0-9.+\-eE]+)", time_line)
    if not match:
        raise ValueError(f"{path} is missing a '# time[fs] = ...' header")

    labels = parse_header_labels(header_line)
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(labels):
        raise ValueError(
            f"{path} has {data.shape[1]} data columns, but header has {len(labels)} labels"
        )
    return float(match.group(1)), labels, data


def label_without_units(label: str) -> str:
    return label.split("[", 1)[0]


def column_index(labels: list[str], column_name: str) -> int:
    lookup = {label: i for i, label in enumerate(labels)}
    if column_name not in lookup:
        bare_lookup = {label_without_units(label): i for i, label in enumerate(labels)}
        if column_name not in bare_lookup:
            choices = ", ".join(labels[1:])
            raise ValueError(
                f"Column {column_name!r} not found. Available data columns: {choices}"
            )
        col = bare_lookup[column_name]
    else:
        col = lookup[column_name]
    if col == 0:
        raise ValueError("The x-face coordinate column cannot be plotted as current")
    return col


def selected_columns(labels: list[str]) -> list[int]:
    requested = getattr(config, "FACE_CURRENT_COLUMNS", None)
    if requested is None:
        return list(range(1, len(labels)))
    if isinstance(requested, str):
        requested = [requested]
    selected = [column_index(labels, name) for name in requested]
    if not selected:
        raise ValueError("FACE_CURRENT_COLUMNS is empty")
    return selected


def load_snapshots(files: list[Path]) -> tuple[list[str], list[int], np.ndarray, np.ndarray, np.ndarray]:
    reference_labels: list[str] | None = None
    reference_x: np.ndarray | None = None
    times: list[float] = []
    currents: list[np.ndarray] = []
    column_indices: list[int] = []

    for path in files:
        time_fs, labels, data = read_face_current(path)
        cols = selected_columns(labels)
        validate_selected_data(path, labels, data, cols)

        if reference_labels is None:
            reference_labels = labels
            reference_x = data[:, 0]
            column_indices = cols
        else:
            if labels != reference_labels:
                raise ValueError(f"{path} uses different columns from {files[0]}")
            if data.shape[0] != reference_x.shape[0] or not np.allclose(data[:, 0], reference_x):
                raise ValueError(f"{path} uses a different x-face grid from {files[0]}")
            if cols != column_indices:
                raise ValueError(f"{path} resolves different current columns from {files[0]}")

        times.append(time_fs)
        currents.append(data[:, cols])

    if reference_labels is None or reference_x is None:
        raise ValueError("No face-current snapshots were loaded")

    order = np.argsort(np.asarray(times))
    sorted_times = np.asarray(times, dtype=float)[order]
    sorted_currents = np.asarray(currents, dtype=float)[order]
    return reference_labels, column_indices, reference_x, sorted_times, sorted_currents


def plot_profiles(
    files: list[Path],
    labels: list[str],
    column_indices: list[int],
    x: np.ndarray,
    times: np.ndarray,
    current: np.ndarray,
    x_limits: tuple[float, float] | None,
) -> None:
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    single_snapshot = len(times) == 1
    for snapshot_idx, time_fs in enumerate(times):
        for series_idx, col in enumerate(column_indices):
            label = labels[col] if single_snapshot else f"{labels[col]}, {time_fs:.2f} fs"
            ax.plot(
                x,
                current[snapshot_idx, :, series_idx],
                linewidth=1.5,
                label=label,
            )

    ax.set_xlabel(labels[0])
    ax.set_ylabel("current density [A/m2]")
    ax.set_title("Face-centered current profiles")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small")

    columns_name = "_".join(filename_token(labels[col]) for col in column_indices)
    output = config.RESULTS_DIR / (
        f"{files[0].stem}_to_{files[-1].stem}_profiles_"
        f"{columns_name}.png"
    )
    save_figure(fig, output, config.DPI)


def main() -> None:
    files = configured_files()
    labels, column_indices, x, times, current = load_snapshots(files)

    x_mask, x_limits = x_range_mask(
        x,
        config.FACE_CURRENT_X_AXIS_RANGE,
        "FACE_CURRENT_X_AXIS_RANGE",
    )
    plot_x = x[x_mask]
    plot_times = times
    plot_current = current[:, x_mask]

    if not np.isfinite(plot_current).all():
        bad_count = int((~np.isfinite(plot_current)).sum())
        raise ValueError(f"Selected face-current data contains {bad_count} non-finite values")

    plot_profiles(files, labels, column_indices, plot_x, plot_times, plot_current, x_limits)


if __name__ == "__main__":
    main()
