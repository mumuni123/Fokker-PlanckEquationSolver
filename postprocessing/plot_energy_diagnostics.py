#!/usr/bin/env python3
"""Plot dKE_bkg and W_bkg_E energy diagnostics vs time from step_diagnostics.dat."""

from __future__ import annotations

from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    filename_token,
    normalize_columns,
    read_table,
    save_figure,
    validate_selected_data,
    x_range_mask,
    y_axis_label,
)

# Columns to plot from step_diagnostics.dat.
ENERGY_DIAG_COLUMNS = [
    "dKE_bkg[J/m2]",
    "W_bkg_E[J/m2]",
]

# Optional x-axis range restriction; set to None to use the full time range,
# or pass (min, max) where either can be None individually.
ENERGY_DIAG_X_AXIS_RANGE = None


def main() -> None:
    diagnostics_file = config.DATA_DIR / "step_diagnostics.dat"
    labels, data = read_table(diagnostics_file)

    column_indices = normalize_columns(labels, ENERGY_DIAG_COLUMNS)
    validate_selected_data(diagnostics_file, labels, data, column_indices)

    time_col = labels.index("time[fs]")
    time = data[:, time_col]

    mask, x_limits = x_range_mask(time, ENERGY_DIAG_X_AXIS_RANGE, "ENERGY_DIAG_X_AXIS_RANGE")
    plot_time = time[mask]

    fig, ax = plt.subplots(figsize=config.FIGSIZE)

    for col in column_indices:
        plot_values = data[mask, col]
        ax.plot(plot_time, plot_values, linewidth=1.8, label=labels[col])

    ax.set_xlabel("time [fs]")
    ax.set_ylabel(y_axis_label(labels, column_indices))
    ax.set_title("Background energy diagnostics")
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    ax.legend()

    columns_name = "_".join(filename_token(labels[col]) for col in column_indices)
    output = config.RESULTS_DIR / f"step_diagnostics_energy_{columns_name}.png"
    save_figure(fig, output, config.DPI)


if __name__ == "__main__":
    main()
