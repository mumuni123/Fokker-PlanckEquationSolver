#!/usr/bin/env python3
"""Plot one or many electric-field spatial profiles."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    figure_output_path,
    filename_token,
    infer_time_settings,
    normalize_columns,
    normalize_plot_columns,
    read_snapshot_table,
    resolve_input_files,
    save_figure,
    snapshot_time_fs,
    validate_selected_data,
    x_range_mask,
    y_axis_label,
)


def _load_snapshot(path: Path) -> tuple[float, list[str], np.ndarray, list[int]]:
    header_time, labels, data = read_snapshot_table(path)
    columns = normalize_columns(labels, config.SPACE_COLUMNS)
    labels, data = normalize_plot_columns(
        labels,
        data,
        columns,
        config.PARAMETERS_FILE,
    )
    validate_selected_data(path, labels, data, columns)
    if header_time is None:
        dt_snapshot_fs, t_end_fs = infer_time_settings(config.PARAMETERS_FILE)
        time_fs = snapshot_time_fs(path, dt_snapshot_fs, t_end_fs)
    else:
        time_fs = header_time
    return time_fs, labels, data, columns


def _plot_single(
    path: Path,
    time_fs: float,
    labels: list[str],
    data: np.ndarray,
    columns: list[int],
) -> None:
    mask, x_limits = x_range_mask(
        data[:, 0],
        config.SPACE_X_AXIS_RANGE,
        "SPACE_X_AXIS_RANGE",
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
    ax.set_title(f"Electric-field profile at t = {time_fs:.3g} fs")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(columns) > 1:
        ax.legend()
    token = "_".join(filename_token(labels[column]) for column in columns)
    save_figure(
        fig,
        config.RESULTS_DIR / f"{path.stem}_{token}.png",
        config.DPI,
    )
    plt.close(fig)


def _plot_overlay(
    files: list[Path],
    snapshots: list[tuple[float, list[str], np.ndarray, list[int]]],
) -> None:
    reference_labels = snapshots[0][1]
    reference_data = snapshots[0][2]
    reference_columns = snapshots[0][3]
    if len(reference_columns) != 1:
        raise ValueError("EX_FIELD_OVERLAY requires exactly one SPACE_COLUMNS entry")
    mask, x_limits = x_range_mask(
        reference_data[:, 0],
        config.SPACE_X_AXIS_RANGE,
        "SPACE_X_AXIS_RANGE",
    )
    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    colors = plt.cm.coolwarm(np.linspace(0.0, 1.0, len(snapshots)))
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
            linewidth=1.5,
            color=color,
            label=f"t = {time_fs:.3g} fs",
        )
    ax.set_xlabel(reference_labels[0])
    ax.set_ylabel(reference_labels[reference_columns[0]])
    ax.set_title("Electric-field evolution")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(snapshots) <= 12:
        ax.legend(title="snapshot", fontsize="small")
    output = figure_output_path(
        config.RESULTS_DIR,
        files,
        "electric_field_evolution",
        filename_token(reference_labels[reference_columns[0]]),
    )
    save_figure(fig, output, config.DPI)
    plt.close(fig)


def main() -> None:
    files = resolve_input_files(config.EX_FIELD_FILES, "EX_FIELD_FILES")
    # ``fields_*.dat`` also matches the staggered-face diagnostics
    # ``fields_face_*.dat``.  Those files contain Ex_face[V/m], not the
    # cell-centred Ex[V/m] requested by this script, so exclude them before
    # loading the batch.
    face_files = [path for path in files if path.name.startswith("fields_face_")]
    files = [path for path in files if not path.name.startswith("fields_face_")]
    if face_files:
        print(
            f"Skipped {len(face_files)} face-centred fields_face_*.dat file(s); "
            "this script processes cell-centred fields_*.dat snapshots."
        )
    if not files:
        raise FileNotFoundError(
            "EX_FIELD_FILES selected only fields_face_*.dat files; "
            "select cell-centred fields_<index>.dat files instead"
        )
    snapshots = [_load_snapshot(path) for path in files]
    print(f"Processing {len(files)} electric-field snapshot(s)")
    for path, snapshot in zip(files, snapshots):
        _plot_single(path, *snapshot)
    if config.EX_FIELD_OVERLAY and len(files) > 1:
        _plot_overlay(files, snapshots)


if __name__ == "__main__":
    main()
