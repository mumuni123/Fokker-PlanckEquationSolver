#!/usr/bin/env python3
"""Plot background electron velocity distribution F(v) from a single snapshot."""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    apply_x_axis_range,
    evaluate_cpp_double_constants,
    filename_token,
    infer_time_settings,
    normalize_columns,
    read_table,
    save_figure,
    snapshot_time_fs,
    validate_selected_data,
    x_range_mask,
)


def light_speed_m_per_s() -> float:
    return evaluate_cpp_double_constants(config.PARAMETERS_FILE).get("c", 2.99792e8)


def u_to_speed_m_per_s(u: np.ndarray, c: float) -> np.ndarray:
    gamma = np.sqrt(1.0 + u * u)
    return c * u / gamma


def fu_to_fv(u: np.ndarray, fu: np.ndarray, c: float) -> np.ndarray:
    gamma = np.sqrt(1.0 + u * u)
    return fu * gamma * gamma * gamma / c


def main() -> None:
    labels, data = read_table(config.FV_FILE)
    column_indices = normalize_columns(labels, config.FV_COLUMNS)
    validate_selected_data(config.FV_FILE, labels, data, column_indices)
    c = light_speed_m_per_s()
    v_m_per_s = u_to_speed_m_per_s(data[:, 0], c)
    plot_columns = {
        col: fu_to_fv(data[:, 0], data[:, col], c)
        for col in column_indices
    }
    mask, x_limits = x_range_mask(
        v_m_per_s,
        config.FV_X_AXIS_RANGE,
        "FV_X_AXIS_RANGE",
    )

    dt_snapshot_fs, t_end_fs = infer_time_settings(config.PARAMETERS_FILE)
    time_fs = snapshot_time_fs(config.FV_FILE, dt_snapshot_fs, t_end_fs)

    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    for col in column_indices:
        ax.plot(v_m_per_s[mask], plot_columns[col][mask], linewidth=1.8, label="F(v)")

    ax.set_xlabel("v[m/s]")
    ax.set_ylabel("F(v)")
    ax.set_title(f"Background electron velocity distribution at t = {time_fs:.2f} fs")
    apply_x_axis_range(ax, x_limits)
    ax.grid(True, alpha=0.3)
    if len(column_indices) > 1:
        ax.legend()

    columns_name = "_".join(filename_token("F(v)") for _ in column_indices)
    output = config.RESULTS_DIR / f"{config.FV_FILE.stem}_{columns_name}.png"
    save_figure(fig, output, config.DPI)


if __name__ == "__main__":
    main()
