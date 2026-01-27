#!/usr/bin/env python3

import argparse
import sys
from glob import glob
from typing import Dict, List, Tuple

import numpy as np

from scipy.optimize import curve_fit
from scipy.special import beta as B

WINDOW_SIZE = 250


def fit_model(x: float | np.ndarray, a: float, b: float) -> float | np.ndarray:
    """
    Fits model double-beta distribution density over ``x`` in [-1, 1].

    The input ``x`` is first mapped from [-1, 1] to [0, 1]. Two asymmetric beta
    components are then combined and scaled to match the expected magnitude for
    the configured ``WINDOW_SIZE``.

    Parameters
    ----------
    x : float | np.ndarray
        Input value(s) in the range [-1, 1].
    a : float
        Fitting parameter for the right component.
    b : float
        Fitting parameter for the left component.

    Returns
    -------
    numpy.ndarray or float
        Scaled density values evaluated at ``x``.
    """
    x_scaled = (x + 1.0) / 2.0

    out = (x_scaled * (1 - x_scaled)**(1 + b)) / B(2, 2 + b)
    out += ((x_scaled)**(1 + a)) / B(2 + a, 2)
    out -= ((x_scaled)**(2 + a)) / B(2 + a, 2)
    amp = (WINDOW_SIZE * 0.1) / 4.0  # Normalization Factor
    out *= amp

    if isinstance(x, np.ndarray):
        for i in range(len(x)):
            if x[i] < -1 or x[i] > 1:
                out[i] = 0
    else:
        if x < -1 or x > 1:
            out = 0

    return out


def load_file(name: str) -> List[Tuple[str, float, float]]:
    """
    Load a single .dat file.

    Parameters
    ----------
    name : str
        Path to the .dat file.

    Returns
    -------
    list
        List of tuples ``(gene, cov_norm, mu)`` for each record.
    """
    out: List[Tuple[str, float, float]] = []

    with open(name, 'r') as f:
        for line in f:
            if line.startswith("#") or "None" in line:
                continue

            line = line.strip().split("\t")
            cov_norm = float(line[0])
            mu = float(line[1])
            # length = int(line[2])
            gene = line[3]
            # strand = line[4].startswith("T")

            out.append((gene, cov_norm, mu))

    return out


def find_state_change(dataset: List[Tuple[str, float, float]]) -> int:
    """
    Finds the first window where the fit indicates a state change.

    Parameters
    ----------
    dataset : list
        Sequence of tuples where the third element is ``mu``.

    Returns
    -------
    int
        Start index of detected change point, or -1 if none found.
    """

    for start in range(len(dataset) // 10):
        start *= 10

        a = [i[2] for i in dataset[start:start+WINDOW_SIZE]]
        bin_count = 20

        bins = np.histogram(a, bins=bin_count, range=(-1, 1))
        x = bins[1][:-1] + (1.0/bin_count)
        bins = bins[0]

        res = run_model(bins, x)
        if res['a'] > 2 or res['b'] > 2:
            print(f"Found! {start} ({res['a']}, {res['b']})", file=sys.stderr)
            return start

    return -1


def sub_file(dataset: List[Tuple[str, float, float]]) -> List[float]:
    """Extract the ``mu`` values from a loaded dataset."""
    return list((i[2] for i in dataset))


def run_model(bins: np.ndarray, x: np.ndarray) -> Dict[str, float]:
    """
    Fit the model ``fit_model`` to histogram ``bins`` at positions ``x``.

    Parameters
    ----------
    bins
        Histogram counts.
    x
        Bin center positions in [-1, 1].

    Returns
    -------
    dict
        Keys ``'a'`` and ``'b'`` with fitted parameters.
    """

    popt, pcov = curve_fit(fit_model, x, bins, bounds=(0, [20, 20]))
    return {'a': popt[0], 'b': popt[1]}


def main(data_path: str) -> None:

    data = load_file(data_path)  # load data
    change = find_state_change(data)  # find state change point

    # if no change was found, use the full dataset
    if change == -1:
        change = len(data)

    data = [i[0] for i in data]
    print("\n".join(data[:change]))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Fit models to find state change points."
    )

    parser.add_argument(
        "dat",
        type=str,
        help="Path to the .dat file to process."
    )

    args = parser.parse_args()
    main(args.dat)
