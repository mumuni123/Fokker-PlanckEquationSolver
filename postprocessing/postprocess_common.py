from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        first = handle.readline().strip()
    if not first.startswith("#"):
        raise ValueError(f"{path} does not start with a header line")
    return first.lstrip("#").split()


def read_table(path: Path) -> tuple[list[str], np.ndarray]:
    labels = read_header(path)
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != len(labels):
        raise ValueError(
            f"{path} has {data.shape[1]} data columns, but header has {len(labels)} labels"
        )
    return labels, data


def evaluate_cpp_double_constants(parameters_file: Path) -> dict[str, float]:
    text = parameters_file.read_text(encoding="utf-8")
    text = re.sub(r"//.*", "", text)
    constants: dict[str, float] = {
        "sqrt": np.sqrt,
        "pow": np.power,
        "pi": np.pi,
    }

    assignment_pattern = re.compile(r"const\s+(?:double|int)\s+(\w+)\s*=\s*(.*?);", re.S)
    for name, expression in assignment_pattern.findall(text):
        py_expression = expression
        py_expression = py_expression.replace("std::sqrt", "sqrt")
        py_expression = py_expression.replace("std::pow", "pow")
        py_expression = re.sub(r"Const::(\w+)", r"\1", py_expression)
        py_expression = re.sub(r"Param::(\w+)", r"\1", py_expression)
        try:
            value = float(eval(py_expression, {"__builtins__": {}}, constants))
        except Exception:
            continue
        constants[name] = value

    return constants


def electric_field_scale(parameters_file: Path) -> float:
    constants = evaluate_cpp_double_constants(parameters_file)
    required = ("qe", "densb", "c", "eps0", "dens", "me")
    missing = [name for name in required if name not in constants]
    if missing:
        raise ValueError(
            f"Cannot compute E0 because {parameters_file} is missing parsable constants: "
            f"{', '.join(missing)}"
        )

    qe = constants["qe"]
    nb = constants["densb"]
    c = constants["c"]
    eps0 = constants["eps0"]
    ne = constants["dens"]
    me = constants["me"]
    omega_p = np.sqrt(ne * qe * qe / (me * eps0))
    return qe * nb * c / (eps0 * omega_p)


def filename_token(text: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.+-]+", "_", text)
    return token.strip("_") or "value"


def normalize_columns(labels: list[str], requested: list[str] | None) -> list[int]:
    if requested is None:
        return list(range(1, len(labels)))

    lookup = {label: i for i, label in enumerate(labels)}
    selected = []
    for name in requested:
        if name not in lookup:
            choices = ", ".join(labels[1:])
            raise ValueError(f"Column {name!r} not found. Available data columns: {choices}")
        if lookup[name] == 0:
            raise ValueError("The x coordinate column cannot be plotted as a dependent variable")
        selected.append(lookup[name])
    return selected


def normalize_plot_columns(
    labels: list[str],
    data: np.ndarray,
    column_indices: list[int],
    parameters_file: Path,
) -> tuple[list[str], np.ndarray]:
    plot_labels = labels.copy()
    plot_data = data.copy()
    e0: float | None = None

    for col in column_indices:
        if labels[col] == "Ex[V/m]":
            if e0 is None:
                e0 = electric_field_scale(parameters_file)
            plot_data[:, col] /= e0
            plot_labels[col] = "Ex/E0"

    return plot_labels, plot_data


def y_axis_label(labels: list[str], column_indices: list[int]) -> str:
    selected = [labels[col] for col in column_indices]
    return ", ".join(selected)


def validate_selected_data(path: Path, labels: list[str], data: np.ndarray, column_indices: list[int]) -> None:
    for col in column_indices:
        values = data[:, col]
        finite = np.isfinite(values)
        if not finite.all():
            bad_count = int((~finite).sum())
            raise ValueError(
                f"{path} column {labels[col]!r} contains {bad_count} non-finite values "
                f"(NaN or inf). Regenerate the simulation output before plotting."
            )


def files_for_prefix(output_dir: Path, prefix: str) -> list[Path]:
    files = sorted(output_dir.glob(f"{prefix}_*.dat"))
    if not files:
        raise FileNotFoundError(f"No files matching {prefix}_*.dat under {output_dir}")
    return files


def parse_first_float_expression(source: str, name: str) -> float | None:
    pattern = re.compile(
        rf"const\s+double\s+{re.escape(name)}\s*=\s*([0-9.+\-eE]+)\s*\*\s*Const::(\w+)"
    )
    match = pattern.search(source)
    if not match:
        return None

    value = float(match.group(1))
    unit = match.group(2)
    scale = {
        "femto": 1.0,
        "micro": 1.0e9,
    }.get(unit)
    if scale is None:
        return None
    return value * scale


def infer_time_settings(parameters_file: Path) -> tuple[float, float | None]:
    if not parameters_file.exists():
        return 0.6, None

    text = parameters_file.read_text(encoding="utf-8")
    dt_snapshot_fs = parse_first_float_expression(text, "dt_snapshot") or 0.6
    t_end_fs = parse_first_float_expression(text, "t_end")
    return dt_snapshot_fs, t_end_fs


def snapshot_index(path: Path) -> int:
    match = re.search(r"_(\d+)\.dat$", path.name)
    if not match:
        raise ValueError(f"Cannot infer snapshot index from {path.name}")
    return int(match.group(1))


def snapshot_time_fs(path: Path, dt_snapshot_fs: float, t_end_fs: float | None) -> float:
    time_fs = snapshot_index(path) * dt_snapshot_fs
    if t_end_fs is not None:
        time_fs = min(time_fs, t_end_fs)
    return time_fs


def interpolate_at_x(
    path: Path,
    x_um: float,
    column_indices: list[int],
    parameters_file: Path,
) -> list[float]:
    labels, data = read_table(path)
    _, data = normalize_plot_columns(labels, data, column_indices, parameters_file)
    validate_selected_data(path, labels, data, column_indices)
    x = data[:, 0]
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if x_um < xmin or x_um > xmax:
        raise ValueError(f"x = {x_um:g} um is outside {path.name} range [{xmin:g}, {xmax:g}] um")
    return [float(np.interp(x_um, x, data[:, col])) for col in column_indices]


def save_figure(fig: plt.Figure, output: Path, dpi: int) -> None:
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    print(f"Saved {output}")
