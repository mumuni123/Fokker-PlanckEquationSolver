#!/usr/bin/env python3
"""Plot beam electron density spatial profile from a single density snapshot."""

from __future__ import annotations

import re

import numpy as np
from matplotlib import pyplot as plt

import config
from postprocess_common import (
    filename_token,
    normalize_columns,
    save_figure,
    validate_selected_data,
    y_axis_label,
)


def read_density_file(path: str):
    """Read a density_*.dat file, returning (time_fs, labels, data)."""
    with open(path, "r", encoding="utf-8") as f:
        time_line = f.readline().strip()
        header_line = f.readline().strip()

    match = re.search(r"time\[fs\]\s*=\s*([0-9.+\-eE]+)", time_line)
    time_fs = float(match.group(1)) if match else 0.0

    labels = header_line.lstrip("#").split()
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return time_fs, labels, data


def main() -> None:
    time_fs, labels, data = read_density_file(str(config.DENSITY_FILE))
    column_indices = normalize_columns(labels, config.DENSITY_COLUMNS)
    validate_selected_data(config.DENSITY_FILE, labels, data, column_indices)

    fig, ax = plt.subplots(figsize=config.FIGSIZE)
    for col in column_indices:
        ax.plot(data[:, 0], data[:, col], linewidth=1.8, label=labels[col])

    ax.set_xlabel(labels[0])
    ax.set_ylabel(y_axis_label(labels, column_indices))
    ax.set_title(f"Beam density profile at t = {time_fs:.2f} fs")
    ax.grid(True, alpha=0.3)
    if len(column_indices) > 1:
        ax.legend()

    columns_name = "_".join(filename_token(labels[col]) for col in column_indices)
    output = config.RESULTS_DIR / f"{config.DENSITY_FILE.stem}_{columns_name}.png"
    save_figure(fig, output, config.DPI)


if __name__ == "__main__":
    main()
