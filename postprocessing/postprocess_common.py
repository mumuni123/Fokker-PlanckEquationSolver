from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_header_labels(header_line: str) -> list[str]:
    labels = re.split(r"\s{2,}", header_line.lstrip("#").strip())
    return labels if len(labels) > 1 else labels[0].split()


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        first = handle.readline().strip()
    if not first.startswith("#"):
        raise ValueError(f"{path} does not start with a header line")
    return parse_header_labels(first)


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


def normalize_x_axis_range(
    range_config: object,
    config_name: str = "X_AXIS_RANGE",
) -> tuple[float | None, float | None] | None:
    if range_config is None:
        return None

    if not isinstance(range_config, (tuple, list)) or len(range_config) != 2:
        raise ValueError(f"{config_name} must be None or a two-item (min, max) range")

    left, right = range_config
    left_value = None if left is None else float(left)
    right_value = None if right is None else float(right)

    if left_value is None and right_value is None:
        return None
    if (
        left_value is not None
        and right_value is not None
        and left_value >= right_value
    ):
        raise ValueError(f"{config_name} lower bound must be smaller than upper bound")

    return left_value, right_value


def x_range_mask(
    x_values: np.ndarray,
    range_config: object,
    config_name: str = "X_AXIS_RANGE",
) -> tuple[np.ndarray, tuple[float, float] | None]:
    x = np.asarray(x_values)
    finite = np.isfinite(x)
    if not finite.all():
        bad_count = int((~finite).sum())
        raise ValueError(f"x axis contains {bad_count} non-finite values")

    normalized = normalize_x_axis_range(range_config, config_name)
    if normalized is None:
        return np.ones(x.shape, dtype=bool), None

    left, right = normalized
    mask = np.ones(x.shape, dtype=bool)
    if left is not None:
        mask &= x >= left
    if right is not None:
        mask &= x <= right

    if not mask.any():
        data_min = float(np.min(x))
        data_max = float(np.max(x))
        requested_left = data_min if left is None else left
        requested_right = data_max if right is None else right
        raise ValueError(
            f"{config_name} range [{requested_left:g}, {requested_right:g}] "
            f"does not overlap data range [{data_min:g}, {data_max:g}]"
        )

    selected_x = x[mask]
    x_limits = (
        float(np.min(selected_x)) if left is None else left,
        float(np.max(selected_x)) if right is None else right,
    )
    return mask, x_limits


def apply_x_axis_range(
    ax: plt.Axes,
    x_limits: tuple[float, float] | None,
) -> None:
    if x_limits is not None:
        ax.set_xlim(left=x_limits[0], right=x_limits[1])


def figure_output_path(
    results_dir: Path,
    files: list[Path],
    plot_name: str,
    columns_name: str,
) -> Path:
    if len(files) == 1:
        file_token = files[0].stem
    else:
        file_token = f"{files[0].stem}_to_{files[-1].stem}"
    return results_dir / f"{file_token}_{plot_name}_{columns_name}.png"


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
