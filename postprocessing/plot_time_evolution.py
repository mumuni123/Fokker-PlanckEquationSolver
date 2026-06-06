#!/usr/bin/env python3
"""Plot full-time evolution at the configured spatial position."""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    filename_token,
    files_for_prefix,
    infer_time_settings,
    interpolate_at_x,
    normalize_plot_columns,
    normalize_columns,
    read_header,
    read_table,
    save_figure,
    snapshot_time_fs,
    validate_selected_data,
    x_range_mask,
    y_axis_label,
)


def main() -> None:
    files = files_for_prefix(config.DATA_DIR, config.TIME_PREFIX)
    labels = read_header(files[0])
    column_indices = normalize_columns(labels, config.TIME_COLUMNS)
    labels, first_data = normalize_plot_columns(
        *read_table(files[0]), column_indices, config.PARAMETERS_FILE
    )
    validate_selected_data(files[0], labels, first_data, column_indices)
    dt_snapshot_fs, t_end_fs = infer_time_settings(config.PARAMETERS_FILE)

    times = np.array([snapshot_time_fs(path, dt_snapshot_fs, t_end_fs) for path in files])
    values = np.array(
        [
            interpolate_at_x(path, config.POSITION_UM, column_indices, config.PARAMETERS_FILE)
            for path in files
        ]
    )
    mask, x_limits = x_range_mask(
        times,
        config.TIME_X_AXIS_RANGE,
        "TIME_X_AXIS_RANGE",
    )
    plot_times = times[mask]
    plot_values = values[mask]
    validation_data = np.column_stack((plot_times, plot_values))
    validation_labels = ["time[fs]", *[labels[col] for col in column_indices]]
    validate_selected_data(
        files[0],
        validation_labels,
        validation_data,
        list(range(1, validation_data.shape[1])),
    )

    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    for series_id, col in enumerate(column_indices):
        ax.plot(plot_times, plot_values[:, series_id], linewidth=1.8, label=labels[col])

    ax.set_xlabel("time [fs]")
    ax.set_ylabel(y_axis_label(labels, column_indices))
    ax.set_title(f"{config.TIME_PREFIX} time evolution at x = {config.POSITION_UM:g} um")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(column_indices) > 1:
        ax.legend()

    columns_name = "_".join(filename_token(labels[col]) for col in column_indices)
    output = (
        config.RESULTS_DIR
        / f"{config.TIME_PREFIX}_time_x{config.POSITION_UM:g}um_{columns_name}.png"
    )
    save_figure(fig, output, config.DPI)


if __name__ == "__main__":
    main()
