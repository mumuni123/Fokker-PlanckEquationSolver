from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

# Simulation output directory that contains density_*.dat and fields_*.dat.
DATA_DIR = ROOT / "output" / "production_120fs_050dt_20260727_173708"

# All postprocessing figures are written here.
RESULTS_DIR = ROOT / "results" / "bkg_density"

# Source parameter file used to infer dt_snapshot, t_end, and electric-field normalization.
PARAMETERS_FILE = ROOT / "src" / "parameters.h"

# Settings for plot_time_evolution.py.
TIME_PREFIX = "fields"
POSITION_UM = 1
TIME_COLUMNS = ["Ex[V/m]"]
TIME_X_AXIS_RANGE = None

# Spatial-profile file selection:
# - one file: DATA_DIR / "fields_00130.dat"
# - explicit batch: [DATA_DIR / "fields_00000.dat", DATA_DIR / "fields_00010.dat"]
# - glob batch: DATA_DIR / "fields_*.dat"
# Every selected snapshot gets its own PNG.  *_OVERLAY additionally creates one
# comparison figure when two or more snapshots are selected.

# Settings for plot_space_profile.py (E_x).
EX_FIELD_FILES = DATA_DIR / "fields_*.dat"
SPACE_COLUMNS = ["Ex[V/m]"]
EX_FIELD_OVERLAY = False
# Use None to plot the full available range, or set (min, max) to restrict it.
# Bounds can also be None individually, e.g. (1.0, None).
SPACE_X_AXIS_RANGE = None

ALL_DENSITY_FILES = "density_*.dat"

# Settings for plot_background_density_evolution.py (n_bkg).
BACKGROUND_DENSITY_FILES = DATA_DIR / ALL_DENSITY_FILES
BACKGROUND_DENSITY_COLUMN = "n_bkg_e[m^-3]"
BACKGROUND_DENSITY_OVERLAY = False
BACKGROUND_DENSITY_X_AXIS_RANGE = None

# Settings for plot_beam_density_profile.py (n_b).
BEAM_DENSITY_FILES = DATA_DIR / ALL_DENSITY_FILES
DENSITY_COLUMNS = ["n_beam[m^-3]"]
BEAM_DENSITY_OVERLAY = False
DENSITY_X_AXIS_RANGE = None

# Settings for plot_fv_profile.py.
FV_FILE = DATA_DIR / "fv_bkg_e_00030.dat"
FV_COLUMNS = ["F(u)"]
FV_X_AXIS_RANGE = (0.0, 3.0e8)

# Settings for plot_current_density_profile.py.
FACE_CURRENT_FILES = DATA_DIR / "current_00025.dat"
# Use None to plot every current column in FACE_CURRENT_FILES.
# Column names may omit units, e.g. "J_total_smoothed" matches "J_total_smoothed[A/m2]".
FACE_CURRENT_COLUMNS = [
    "J_total_ampere_face",
]
FACE_CURRENT_X_AXIS_RANGE = None

# Settings for plot_electron_kinetics.py.
# Output is kept in one task-specific directory, like the EPOCH postprocessor.
KINETIC_RESULTS_DIR = ROOT / "results" / DATA_DIR.name / "electron_kinetics"
# The solver stores 1D energies per transverse area.  Multiplication by this
# physical/implicit area converts every plotted integrated energy to joules.
KINETIC_TRANSVERSE_AREA_M2 = 1.0
# Requested comparison times.  The nearest available snapshots are selected;
# repeated nearest files (for times beyond an incomplete run) are removed.
KINETIC_SAMPLE_TIMES_FS = [0.0, 25.0, 50.0, 75.0, 100.0, 120.0]
# Distribution representation: "count", "density", or "probability".
KINETIC_DISTRIBUTION_MODE = "density"
# Explicit file lists override time-based automatic selection.
KINETIC_SPECTRUM_FILES = None
KINETIC_PARALLEL_MOMENTUM_FILES = None
KINETIC_PERPENDICULAR_MOMENTUM_FILES = None
# Empty/zero bins are omitted on logarithmic axes.
KINETIC_LOG_Y = True
# Plotting floor relative to the snapshot-integrated particle number.  The raw
# CSV retains every bin; this only hides floating-point underflow tails that
# otherwise span hundreds of decades on a logarithmic axis.
KINETIC_MIN_BIN_FRACTION = 1.0e-16
# Conservative neighboring-bin averaging for the energy-spectrum PNG only.
# CSV files and all integral summaries always retain the unsmoothed data.
KINETIC_SPECTRUM_SMOOTHING_BINS = 5
# Settings for plot_electron_energy_spectra.py.  The two plots use the same
# conservative global spectrum but expose the thermal/bulk and nonthermal-tail
# windows separately.  ``None`` means the largest available bin edge.
KINETIC_LOW_ENERGY_RANGE_EV = (1.0, 2.0e4)
KINETIC_HIGH_ENERGY_RANGE_EV = (2.0e4, None)

# Plot style and output quality.
FIGSIZE = (7.2, 4.5)
DPI = 200
