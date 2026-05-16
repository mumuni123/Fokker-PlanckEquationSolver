#!/usr/bin/env python3
"""Plot full-space profiles from the configured snapshot file."""

from __future__ import annotations

from matplotlib import pyplot as plt

import config
from postprocess_common import (
    filename_token,
    normalize_plot_columns,
    normalize_columns,
    read_table,
    save_figure,
    validate_selected_data,
    y_axis_label,
)


def main() -> None:
    labels, data = read_table(config.SPACE_FILE)
    column_indices = normalize_columns(labels, config.SPACE_COLUMNS)
    labels, data = normalize_plot_columns(
        labels, data, column_indices, config.PARAMETERS_FILE
    )
    validate_selected_data(config.SPACE_FILE, labels, data, column_indices)

    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    for col in column_indices:
        ax.plot(data[:, 0], data[:, col], linewidth=1.8, label=labels[col])

    ax.set_xlabel(labels[0])
    ax.set_ylabel(y_axis_label(labels, column_indices))
    ax.set_title(f"Spatial profile: {config.SPACE_FILE.name}")
    ax.grid(True, alpha=0.3)
    if len(column_indices) > 1:
        ax.legend()

    columns_name = "_".join(filename_token(labels[col]) for col in column_indices)
    output = config.RESULTS_DIR / f"{config.SPACE_FILE.stem}_{columns_name}.png"
    save_figure(fig, output, config.DPI)


if __name__ == "__main__":
    main()
