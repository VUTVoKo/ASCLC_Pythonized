"""Notebook-facing loading and export functions; no example file dependencies."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from asclc_backend import Measurements, load_measurements, current_density


@dataclass
class Calculation:
    measurements: Measurements
    V: np.ndarray
    J: np.ndarray


def run_calculation(path, *, device):
    """Load the standard CSV and return measured voltage and current density.

    ``V = U + device['voltage_offset']`` and ``J = I / device['area']``.
    Signs, row order and nonfinite values are preserved; ``device`` is not mutated.
    """
    measurements = load_measurements(path)
    if not np.isfinite(device['voltage_offset']):
        raise ValueError("Voltage offset must be finite.")
    v = measurements.U + device['voltage_offset']
    j = current_density(measurements.I, device['area'])
    return Calculation(measurements, v, j)


def export_results(result, figures, directory):
    """Save each named figure as PNG/PDF/SVG plus the measured CSV.

    ``figures`` maps a base filename (no extension) to a Matplotlib figure.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, figure in figures.items():
        for extension in ('png', 'pdf', 'svg'):
            figure.savefig(directory / f'{name}.{extension}', dpi=240,
                           bbox_inches='tight')
    np.savetxt(directory / 'measured.csv',
               np.column_stack((result.measurements.U, result.measurements.I,
                                result.V, result.J)),
               delimiter=',', comments='',
               header='U_V,I_A,V_V,J_Am2')
    return directory
