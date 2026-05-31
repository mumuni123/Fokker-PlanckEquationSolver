#!/usr/bin/env python3
"""Plot background electron velocity distribution F(v) from a single snapshot."""

from __future__ import annotations

from matplotlib import pyplot as plt

import config
from postprocess_common import (
    filename_token,
    infer_time_settings,
    normalize_columns,
    read_table,
    save_figure,
    snapshot_index,
    snapshot_time_fs,
    validate_selected_data,
    y_axis_label,
)


def main() -> None:
    labels, data = read_table(config.FV_FILE)
    column_indices = normalize_columns(labels, config.FV_COLUMNS)
    validate_selected_data(config.FV_FILE, labels, data, column_indices)

    idx = snapshot_index(config.FV_FILE)
    dt_snapshot_fs, t_end_fs = infer_time_settings(config.PARAMETERS_FILE)
    time_fs = snapshot_time_fs(config.FV_FILE, dt_snapshot_fs, t_end_fs)

    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    for col in column_indices:
        ax.plot(data[:, 0], data[:, col], linewidth=1.8, label=labels[col])

    ax.set_xlabel(labels[0])
    ax.set_ylabel(y_axis_label(labels, column_indices))
    ax.set_title(f"Background electron velocity distribution at t = {time_fs:.2f} fs")
    ax.grid(True, alpha=0.3)
    if len(column_indices) > 1:
        ax.legend()

    columns_name = "_".join(filename_token(labels[col]) for col in column_indices)
    output = config.RESULTS_DIR / f"{config.FV_FILE.stem}_{columns_name}.png"
    save_figure(fig, output, config.DPI)


if __name__ == "__main__":
    main()
