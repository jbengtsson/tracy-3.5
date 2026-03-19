from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from enum import IntEnum
from pathlib import Path
from typing import TextIO

import numpy as np
from numpy.typing import NDArray


class Plane(IntEnum):
    X = 0
    Y = 1
    Z = 2


class PhaseSpace(IntEnum):
    X = 0
    PX = 1
    Y = 2
    PY = 3


X_ = Plane.X
Y_ = Plane.Y
Z_ = Plane.Z

x_ = PhaseSpace.X
px_ = PhaseSpace.PX
y_ = PhaseSpace.Y
py_ = PhaseSpace.PY


PlaneArray = NDArray[np.float64]


def sqr(x: float | NDArray[np.float64]) -> float | NDArray[np.float64]:
    return x**2


def print_write(fmt: str, *args: object) -> None:
    print(fmt % args, end="")


def file_write(outf: TextIO, fmt: str, *args: object) -> None:
    outf.write(fmt % args)



@dataclass
class LinOpt:
    loc: list[int] = field(default_factory=list)
    name: list[str] = field(default_factory=list)
    s: PlaneArray = field(default_factory=lambda: np.zeros(0, dtype=float))
    alpha: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    beta: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    nu: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    eta: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    etap: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))

    def read(self, file_name: str | Path, verbose: bool = False) -> None:
        rows: list[tuple[int, str, float, PlaneArray, PlaneArray, PlaneArray,
                         PlaneArray, PlaneArray]] = []

        with open(file_name, "r", encoding="utf-8") as inf:
            for raw_line in inf:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue

                tokens = line.split()
                if len(tokens) < 14:
                    raise ValueError(f"Malformed linlat row: {line}")
                n = int(tokens[0])
                name = tokens[1].strip()
                s = float(tokens[2])

                alpha = \
                    np.array([float(tokens[4]), float(tokens[9])], dtype=float)
                beta = \
                    np.array([float(tokens[5]), float(tokens[10])], dtype=float)
                nu = \
                    np.array([float(tokens[6]), float(tokens[11])], dtype=float)
                eta = \
                    np.array([float(tokens[7]), float(tokens[12])], dtype=float)
                etap = \
                    np.array([float(tokens[8]), float(tokens[13])], dtype=float)
                rows.append((n, name, s, alpha, beta, nu, eta, etap))

                if verbose:
                    print_write(
                        "%4d, %-15s, %9.5f, %9.5f, %8.5f, %8.5f, %8.5f, %8.5f"
                        ", %9.5f, %8.5f, %8.5f, %8.5f, %8.5f\n",
                        n,
                        name,
                        s,
                        alpha[X_],
                        beta[X_],
                        nu[X_],
                        eta[X_],
                        etap[X_],
                        alpha[Y_],
                        beta[Y_],
                        nu[Y_],
                        eta[Y_],
                        etap[Y_],
                    )

        self.loc = [row[0] for row in rows]
        self.name = [row[1] for row in rows]
        self.s = np.array([row[2] for row in rows], dtype=float)
        self.alpha = \
            np.column_stack([row[3] for row in rows]) \
            if rows else np.zeros((2, 0), dtype=float)
        self.beta = \
            np.column_stack([row[4] for row in rows]) \
            if rows else np.zeros((2, 0), dtype=float)
        self.nu = \
            np.column_stack([row[5] for row in rows]) \
            if rows else np.zeros((2, 0), dtype=float)
        self.eta = \
            np.column_stack([row[6] for row in rows]) \
            if rows else np.zeros((2, 0), dtype=float)
        self.etap = \
            np.column_stack([row[7] for row in rows]) \
            if rows else np.zeros((2, 0), dtype=float)


@dataclass
class BPMData:
    n_bpm: int = 0
    n_turn: int = 0
    name: list[str] = field(default_factory=list)
    loc: list[int] = field(default_factory=list)
    data: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0, 0), dtype=float))

    @staticmethod
    def normalize_bpm_name(name: str) -> str:
        return re.sub("-", "_", name).lower()

    def get_loc(self, name: str, lin_opt: LinOpt) -> int:
        try:
            return lin_opt.name.index(name)
        except ValueError as exc:
            raise ValueError(f"BPM '{name}' not found in linear optics table") \
                from exc

    def read_bpm_names(
            self, inf: TextIO, lin_opt: LinOpt, verbose: bool = False) -> None:
        n_print = 8
        inf.readline()
        header = inf.readline().strip().split()
        if len(header) < 2:
            raise ValueError(
                "Malformed BPM header: expected number of BPMs and turns")
        self.n_bpm, self.n_turn = map(int, header[:2])
        inf.readline()
        print_write(
            "\nno of BPMs = %d, no of turns = %d \n", self.n_bpm, self.n_turn)

        self.data = np.zeros((2, self.n_bpm, self.n_turn), dtype=float)
        self.name.clear()
        self.loc.clear()

        if verbose:
            print_write("\n")

        for j in range(self.n_bpm):
            parts = inf.readline().strip().split()
            if not parts:
                raise ValueError("Malformed BPM name row")
            bpm_name = self.normalize_bpm_name(parts[0])
            self.name.append(bpm_name)
            self.loc.append(self.get_loc(bpm_name, lin_opt))
            if verbose:
                print_write(" %s", self.name[j])
                if (j + 1) % n_print == 0:
                    print_write("\n")

        if verbose and (self.n_bpm % n_print != 0):
            print_write("\n")

    def read_bpm_data(
            self, plane: int, inf: TextIO, verbose: bool = False) -> None:
        n_print = 8
        if verbose:
            print_write("\n")

        inf.readline()
        for k in range(self.n_turn):
            if verbose:
                print_write("\n")

            data = inf.readline().strip().split()
            tail = inf.readline().strip().split()
            if len(data) != self.n_bpm - 1 or not tail:
                raise ValueError("Malformed TbT row")
            data.append(tail[0])
            self.data[plane, :, k] = 1e-3 * np.array(data, dtype=float)

            if verbose:
                for j in range(self.n_bpm):
                    print_write("%10.6f", 1e3 * self.data[plane, j, k])
                    if (j + 1) % n_print == 0:
                        print_write("\n")
                if self.n_bpm % n_print != 0:
                    print_write("\n")

    def read_tbt(self, file_name: str | Path, lin_opt: LinOpt) -> None:
        with open(file_name, "r", encoding="utf-8") as inf:
            self.read_bpm_names(inf, lin_opt)
            # assert False
            self.read_bpm_data(X_, inf)
            self.read_bpm_data(Y_, inf)


@dataclass
class EstimatedLinOpt:
    beta_pinger: PlaneArray = \
        field(default_factory=lambda: np.zeros(2, dtype=float))
    alpha_mean: PlaneArray = \
        field(default_factory=lambda: np.zeros(2, dtype=float))
    alpha_sigma: PlaneArray = \
        field(default_factory=lambda: np.zeros(2, dtype=float))
    tune_mean: PlaneArray = \
        field(default_factory=lambda: np.zeros(2, dtype=float))
    tune_sigma: PlaneArray = \
        field(default_factory=lambda: np.zeros(2, dtype=float))
    beta: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))

    beta_mean: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    beta_sigma: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    nu: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))

    dnu_mean: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    dnu_sigma: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    beta_samples: list[PlaneArray] = field(default_factory=list)
    dnu_samples: list[PlaneArray] = field(default_factory=list)
    twoJ: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    phi: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))
    phi0: PlaneArray = \
        field(default_factory=lambda: np.zeros((2, 0), dtype=float))

    def zero(self, n: int) -> None:
        shape = (2, n)
        self.beta = np.zeros(shape, dtype=float)
        self.beta_mean = np.zeros(shape, dtype=float)
        self.beta_sigma = np.zeros(shape, dtype=float)
        self.nu = np.zeros(shape, dtype=float)
        self.dnu_mean = np.zeros(shape, dtype=float)
        self.dnu_sigma = np.zeros(shape, dtype=float)
        self.beta_samples = []
        self.dnu_samples = []

    def get_stats(
            self, bpm_data: BPMData, lin_opt: LinOpt,
            out_file: str | Path = "lin_opt.out") -> None:
        dbeta_max = 5.0
        dnu_max = 0.05

        self.beta_mean = np.zeros((2, bpm_data.n_bpm), dtype=float)
        self.beta_sigma = np.zeros((2, bpm_data.n_bpm), dtype=float)
        self.dnu_mean = np.zeros((2, bpm_data.n_bpm), dtype=float)
        self.dnu_sigma = np.zeros((2, bpm_data.n_bpm), dtype=float)

        if not self.beta_samples or not self.dnu_samples:
            raise ValueError("No statistics samples accumulated")

        beta_samples = np.stack(self.beta_samples, axis=0)
        dnu_samples = np.stack(self.dnu_samples, axis=0)
        self.beta_mean = np.mean(beta_samples, axis=0)
        self.dnu_mean = np.mean(dnu_samples, axis=0)

        if beta_samples.shape[0] > 1:
            self.beta_sigma = np.std(beta_samples, axis=0, ddof=1)
            self.dnu_sigma = np.std(dnu_samples, axis=0, ddof=1)
        else:
            self.beta_sigma = np.zeros_like(self.beta_mean)
            self.dnu_sigma = np.zeros_like(self.dnu_mean)

        with open(out_file, "w", encoding="utf-8") as outf:
            dbeta = np.zeros(2, dtype=float)
            dnu = np.zeros(2, dtype=float)
            file_write(
                outf, "\n# bpm  s [m]                 beta [m]"
                "                           nu\n")
            for j in range(bpm_data.n_bpm):
                loc = bpm_data.loc[j]
                for k in range(2):
                    dbeta[k] = self.beta_mean[k, j]
                    if self.beta_sigma[k, j] > dbeta_max:
                        dbeta[k] = 0.0
                        self.beta_sigma[k, j] = 0.0

                    dnu[k] = \
                        self.dnu_mean[k, j] \
                        - (lin_opt.nu[k, loc] - int(lin_opt.nu[k, loc]))
                    if self.dnu_sigma[k, j] > dnu_max:
                        print_write("\nBPM # %d excluded, plane = %d\n", j, k)
                        dnu[k] = 0.0
                        self.dnu_sigma[k, j] = 0.0

                file_write(
                    outf,
                    "%4d %8.3f %7.3f +/- %5.3f %7.3f +/- %5.3f%7.3f +/- %5.3f"
                    " %7.3f +/- %5.3f %8.3f %8.3f\n",
                    j + 1,
                    lin_opt.s[loc],
                    dbeta[X_],
                    self.beta_sigma[X_, j],
                    dbeta[Y_],
                    self.beta_sigma[Y_, j],
                    self.dnu_mean[X_, j],
                    self.dnu_sigma[X_, j],
                    self.dnu_mean[Y_, j],
                    self.dnu_sigma[Y_, j],
                    lin_opt.beta[X_, loc],
                    lin_opt.beta[Y_, loc],
                )


def get_window_weights(n: int, window: int) -> NDArray[np.float64]:
    """Return the selected analysis window of length n."""
    if window == 1:
        return np.ones(n, dtype=float)
    if window == 2:
        return np.sin(np.arange(n, dtype=float) / float(n - 1) * np.pi)
    if window == 3:
        return np.sin(np.arange(n, dtype=float) / float(n - 1) * np.pi) ** 2
    raise ValueError(f"Unsupported window type: {window}")


def coherent_gain(window: int, n: int) -> float:
    """Return the coherent gain for the selected window."""
    return float(np.mean(get_window_weights(n, window)))


def apply_window(x: NDArray[np.float64], window: int) -> NDArray[np.float64]:
    return x * get_window_weights(x.shape[-1], window)


def FFT1(x: PlaneArray, window: int) -> tuple[PlaneArray, PlaneArray]:
    """Return amplitude and phase spectra with coherent-gain correction."""
    n = len(x)
    x1 = np.fft.rfft(apply_window(x, window))
    gain = coherent_gain(window, n)
    return np.abs(x1) * 2.0 / (n * gain), np.angle(x1)


def FFT2(x: PlaneArray, window: int) -> NDArray[np.complex128]:
    """Return the complex FFT spectrum for the selected window."""
    return np.fft.rfft(apply_window(x, window))


def get_ind(n: int, k: int) -> tuple[int, int]:
    if k == 0:
        return 1, 1
    if k == n // 2:
        return n // 2 - 1, n // 2 - 1
    return k - 1, k + 1


def get_nu1(n: int, A: PlaneArray, k: int, window: int) -> float:
    """Estimate the tune from the peak bin and its dominant neighbor."""
    ind1, ind3 = get_ind(n, k)
    if A[ind3] > A[ind1]:
        A1, A2, ind = A[k], A[ind3], k
    else:
        A1, A2 = A[ind1], A[k]
        ind = ind1 if k != 0 else -1

    if A1 + A2 == 0.0:
        return 0.0
    if window == 1:
        return (ind + A2 / (A1 + A2)) / n
    if window == 2:
        return (ind - 0.5 + 2.0 * A2 / (A1 + A2)) / n
    if window == 3:
        return (ind - 1.0 + 3.0 * A2 / (A1 + A2)) / n
    raise ValueError(f"Unsupported window type: {window}")


def sinc(omega: float) -> float:
    return math.sin(omega) / omega if omega != 0.0 else 1.0


def get_A(n: int, A: PlaneArray, nu: float, k: int, window: int) -> float:
    """Recover the line amplitude from the interpolated tune."""
    if window == 1:
        corr = sinc(np.pi * (k - nu * n))
    elif window == 2:
        corr = \
            (sinc(np.pi * (k + 0.5 - nu * n)) \
             + sinc(np.pi * (k - 0.5 - nu * n))) / 2.0
    elif window == 3:
        raise NotImplementedError("get_A is not implemented for window=3")
    else:
        raise ValueError(f"Unsupported window type: {window}")
    if abs(corr) < 1e-14:
        return 0.0
    return A[k] / corr


def get_alpha(
        n: int, X: NDArray[np.complex128],
        nu: float, k: int) -> tuple[float, float]:
    """Estimate the local spectral phase shift and damping parameter."""
    I = complex(0.0, 1.0)
    ind1, ind3 = get_ind(n, k)
    if abs(X[ind3]) > abs(X[ind1]):
        d, rho = 1, X[ind3] / X[k]
    else:
        d, rho = -1, X[ind1] / X[k]
    z = \
        (1.0 - rho) \
        / (1.0 - rho * np.exp(-I * 2.0 * np.pi * float(d) / float(n)))
    delta = n * np.angle(z) / (2.0 * np.pi)
    alpha = n * math.log(abs(z)) / (2.0 * np.pi)
    return delta, alpha


def get_peak(n: int, A: PlaneArray) -> int:
    """Return the dominant local-maximum FFT bin."""
    k = 0
    peak = 0.0
    for ind2 in range(n // 2 + 1):
        ind1, ind3 = get_ind(n, ind2)
        if A[ind2] > peak and A[ind1] < A[ind2] and A[ind2] > A[ind3]:
            peak = A[ind2]
            k = ind2
    return k


def get_phi(n: int, k: int, nu: float, phi: PlaneArray) -> float:
    """Return the phase corrected to the interpolated tune."""
    phi_nu = phi[k] - (n * nu - k) * np.pi
    if phi_nu > np.pi:
        phi_nu -= 2.0 * np.pi
    elif phi_nu < -np.pi:
        phi_nu += 2.0 * np.pi
    return phi_nu


def get_nu2(
        n: int, x: PlaneArray,
        window: int) -> tuple[float, float, float, float, float]:
    """Estimate tune, amplitude, phase, phase shift"
    ", and damping from TbT data."""
    A, phi = FFT1(x, window)
    x_fft = FFT2(x, 1)
    return get_nu2_from_spectra(n, A, phi, x_fft, window)


def get_nu2_from_spectra(
    n: int,
    A: PlaneArray,
    phi: PlaneArray,
    x_fft: NDArray[np.complex128],
    window: int,
) -> tuple[float, float, float, float, float]:
    """Apply the estimator chain to precomputed spectra."""
    k = get_peak(n, A)
    nu = get_nu1(n, A, k, window)
    A_nu = get_A(n, A, nu, k, window)
    phi_nu = get_phi(n, k, nu, phi)
    delta, alpha = get_alpha(n, x_fft, nu, k)
    return nu, A_nu, phi_nu, delta, alpha



def get_nus(
    outf: TextIO,
    cut: int,
    n: int,
    window: int,
    bpm_data: BPMData,
    lin_opt: LinOpt,
    est_lin_opt: EstimatedLinOpt,
) -> None:
    sgn = [1, -1]

    tune_sum = np.zeros(2, dtype=float)
    tune_sum2 = np.zeros(2, dtype=float)
    alpha_sum = np.zeros(2, dtype=float)
    alpha_sum2 = np.zeros(2, dtype=float)
    twoJ_sum = np.zeros(2, dtype=float)
    twoJ_sum2 = np.zeros(2, dtype=float)
    phi0_sum = np.zeros(2, dtype=float)
    phi0_sum2 = np.zeros(2, dtype=float)

    tunes = np.zeros((bpm_data.n_bpm, 2), dtype=float)
    As = np.zeros((bpm_data.n_bpm, 2), dtype=float)
    phis = np.zeros((bpm_data.n_bpm, 2), dtype=float)
    nus = np.zeros((bpm_data.n_bpm, 2), dtype=float)
    delta = np.zeros(2, dtype=float)
    alpha = np.zeros(2, dtype=float)
    phi0 = np.zeros(2, dtype=float)

    segment = bpm_data.data[:, :, cut : n + cut].copy()
    segment -= np.mean(segment, axis=2, keepdims=True)

    spectra_amp = np.zeros((2, bpm_data.n_bpm, n // 2 + 1), dtype=float)
    spectra_phi = np.zeros((2, bpm_data.n_bpm, n // 2 + 1), dtype=float)
    spectra_rect = \
        np.zeros((2, bpm_data.n_bpm, n // 2 + 1), dtype=np.complex128)

    for j in range(2):
        fft_windowed = np.fft.rfft(apply_window(segment[j], window), axis=1)
        gain = coherent_gain(window, n)
        spectra_amp[j] = np.abs(fft_windowed) * 2.0 / (n * gain)
        spectra_phi[j] = np.angle(fft_windowed)
                # Rectangular FFT kept for alpha estimation.
        spectra_rect[j] = np.fft.rfft(apply_window(segment[j], 1), axis=1)

    print_write("\n")
    for i in range(bpm_data.n_bpm):
        loc = bpm_data.loc[i]
        for j in range(2):
            tunes[i, j], As[i, j], phis[i, j], delta[j], alpha[j] = \
                get_nu2_from_spectra(
                n,
                spectra_amp[j, i],
                spectra_phi[j, i],
                spectra_rect[j, i],
                window,
            )

            if sgn[j] < 0:
                phis[i, j] = -phis[i, j]
            if phis[i, j] < 0.0:
                phis[i, j] += 2.0 * np.pi
            nus[i, j] = phis[i, j] / (2.0 * np.pi)

            tune_sum[j] += tunes[i, j]
            tune_sum2[j] += sqr(tunes[i, j])
            alpha_sum[j] += alpha[j]
            alpha_sum2[j] += sqr(alpha[j])

            twoJ = sqr(As[i, j]) / lin_opt.beta[j, loc]
            twoJ_sum[j] += twoJ
            twoJ_sum2[j] += sqr(twoJ)

            phi0[j] = (nus[i, j] - (lin_opt.nu[j, loc] % 1.0)) * 2.0 * np.pi
            if phi0[j] < 0.0:
                phi0[j] += 2.0 * np.pi
            phi0_sum[j] += phi0[j]
            phi0_sum2[j] += sqr(phi0[j])

    twoJ_mean = twoJ_sum / bpm_data.n_bpm
    twoJ_sigma = \
        np.sqrt((bpm_data.n_bpm * twoJ_sum2 - sqr(twoJ_sum)) / (bpm_data.n_bpm * (bpm_data.n_bpm - 1.0)))

    phi0_mean = phi0_sum / bpm_data.n_bpm
    phi0_sigma = \
        np.sqrt((bpm_data.n_bpm * phi0_sum2 - np.square(phi0_sum)) \
                / (bpm_data.n_bpm * (bpm_data.n_bpm - 1.0)))

    print_write(
        "\ntwoJ  = [%9.3e+/-%9.3e, %9.3e+/-%9.3e]"
        ", phi0 = [%5.3f+/-%5.3f, %5.3f+/-%5.3f]\n",
        twoJ_mean[X_],
        twoJ_sigma[X_],
        twoJ_mean[Y_],
        twoJ_sigma[Y_],
        phi0_mean[X_],
        phi0_sigma[X_],
        phi0_mean[Y_],
        phi0_sigma[Y_],
    )
    print_write(
        "A0    = [%5.3f, %5.3f] mm\n",
        1e3 * math.sqrt(twoJ_mean[X_] * est_lin_opt.beta_pinger[X_]),
        1e3 * math.sqrt(twoJ_mean[Y_] * est_lin_opt.beta_pinger[Y_]),
    )

    dnu = np.zeros(2, dtype=float)
    beta_run = np.zeros((2, bpm_data.n_bpm), dtype=float)
    dnu_run = np.zeros((2, bpm_data.n_bpm), dtype=float)
    for i in range(bpm_data.n_bpm):
        loc = bpm_data.loc[i]
        for j in range(2):
            beta = sqr(As[i, j]) / twoJ_mean[j]
            nus[i, j] -= phi0_mean[j] / (2.0 * np.pi)
            if nus[i, j] < 0.0:
                nus[i, j] += 1.0

            dnu[j] = nus[i, j] - (lin_opt.nu[j, loc] % 1.0)
            if dnu[j] < -0.5:
                dnu[j] += 1.0
            if dnu[j] > 0.5:
                dnu[j] -= 1.0

            beta_run[j, i] = beta
            dnu_run[j, i] = dnu[j]

        file_write(
            outf,
            "%4d %7.3f %8.3f %8.3f\n", i + 1, lin_opt.s[loc], dnu[X_], dnu[Y_])

    est_lin_opt.beta_samples.append(beta_run)
    est_lin_opt.dnu_samples.append(dnu_run)

    for j in range(2):
        est_lin_opt.tune_mean[j] = tune_sum[j] / bpm_data.n_bpm
        if sgn[j] < 0:
            est_lin_opt.tune_mean[j] = 1.0 - est_lin_opt.tune_mean[j]
        est_lin_opt.tune_sigma[j] = math.sqrt(
            (bpm_data.n_bpm * tune_sum2[j] - sqr(tune_sum[j])) \
            / (bpm_data.n_bpm * (bpm_data.n_bpm - 1.0))
        )

        est_lin_opt.alpha_mean[j] = alpha_sum[j] / bpm_data.n_bpm
        est_lin_opt.alpha_sigma[j] = math.sqrt(
            (bpm_data.n_bpm * alpha_sum2[j] - sqr(alpha_sum[j])) \
            / (bpm_data.n_bpm * (bpm_data.n_bpm - 1.0))
        )

    print_write(
        "\nnu    = [%9.6f+/-%8.6f, %9.6f+/-%8.6f]\n",
        est_lin_opt.tune_mean[X_],
        est_lin_opt.tune_sigma[X_],
        est_lin_opt.tune_mean[Y_],
        est_lin_opt.tune_sigma[Y_],
    )
    print_write(
        "alpha = [%9.6f+/-%8.6f, %9.6f+/-%8.6f]\n",
        est_lin_opt.alpha_mean[X_],
        est_lin_opt.alpha_sigma[X_],
        est_lin_opt.alpha_mean[Y_],
        est_lin_opt.alpha_sigma[Y_],
    )
    print_write(
        "%8.5f %8.5f\n", nus[6, X_] - nus[5, X_], nus[6, Y_] - nus[5, Y_])


def prt_FFT(cut: int, xy: PlaneArray, window: int) -> None:
    n = len(xy[X_])
    with open("tbt_data.out", "w", encoding="utf-8") as outf:
        for j in range(cut, n + cut):
            file_write(outf, "%5d %11.3e %11.3e\n", j + 1, xy[X_, j], xy[Y_, j])

    x1 = xy[:, cut : n + cut]
    A = np.zeros((2, n // 2 + 1), dtype=float)
    phi = np.zeros((2, n // 2 + 1), dtype=float)
    for k in range(2):
        A[k], phi[k] = FFT1(x1[k], window)

    with open("tbt_fft.out", "w", encoding="utf-8") as outf:
        for k in range(n // 2 + 1):
            file_write(
                outf, "%5d %9.3e %9.3e %9.3e\n",
                k + 1, float(k) / float(n), A[X_, k], A[Y_, k])


def get_b1ob2_dnu(
        n: int, ps1: PlaneArray,
        ps2: PlaneArray) -> tuple[PlaneArray, PlaneArray]:
    """Estimate beta ratio and phase advance from two phase-space traces."""
    print_write("\n")
    b1ob2 = np.zeros(2, dtype=float)
    dnu = np.zeros(2, dtype=float)
    for k in range(2):
        x1_sqr = np.sum(sqr(ps1[k])) / n
        x2_sqr = np.sum(sqr(ps2[k])) / n
        x1x2 = np.sum(ps1[k] * ps2[k]) / n
        b1ob2[k] = x1_sqr / x2_sqr
        dnu[k] = math.acos(x1x2 / math.sqrt(x1_sqr * x2_sqr)) / (2.0 * np.pi)
        print_write("b1ob2 = %9.3e, dnu = %5.3f\n", b1ob2[k], dnu[k])
    return b1ob2, dnu


def ss_est(
        cut: int, n: int, bpm1: int, bpm2: int,
        bpm_data: BPMData, lin_opt: LinOpt) -> None:
    ps1 = bpm_data.data[:, bpm1 - 1, cut : n + cut].copy()
    ps2 = bpm_data.data[:, bpm2 - 1, cut : n + cut].copy()

    get_b1ob2_dnu(n, ps1, ps2)

    loc1 = bpm_data.loc[bpm1 - 1]
    loc2 = bpm_data.loc[bpm2 - 1]
    beta1 = lin_opt.beta[:, loc1]
    b1ob2 = lin_opt.beta[:, loc1] / lin_opt.beta[:, loc2]
    dnu = lin_opt.nu[:, loc2] - lin_opt.nu[:, loc1]
    print_write("\n")
    for k in range(2):
        print_write("b1ob2 = %9.3e, dnu = %5.3f\n", b1ob2[k], dnu[k])

    with open("tbt_phase_space.out", "w", encoding="utf-8") as outf:
        ps = np.zeros(4, dtype=float)
        twoJ = np.zeros(2, dtype=float)
        for j in range(n):
            for k in range(2):
                ps[2 * k] = ps1[k, j] / math.sqrt(beta1[k])
                ps[2 * k + 1] = (
                    math.sqrt(b1ob2[k]) * ps2[k, j] - ps1[k, j] \
                    * math.cos(2.0 * math.pi * dnu[k])
                ) / (math.sqrt(beta1[k]) * math.sin(2.0 * math.pi * dnu[k]))
                twoJ[k] = sqr(ps[2 * k]) + sqr(ps[2 * k + 1])
            file_write(
                outf,
                "%4d%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n",
                j + 1,
                ps[x_],
                ps[px_],
                ps[y_],
                ps[py_],
                twoJ[X_],
                twoJ[Y_],
            )


def prt_name(outf: TextIO, name: str) -> None:
    parts = name.split(" ", 1)
    file_write(outf, "%s," % parts[0])
    if len(parts) > 1:
        file_write(outf, "%s" % parts[1])


def run(home_dir: str | Path) -> None:
    home = Path(home_dir)
    bpm_data = BPMData()
    lin_opt = LinOpt()
    est_lin_opt = EstimatedLinOpt()

    lin_opt.read(home / "linlat.out")

    window = 2
    cut = 0
    n_turn = 2 * 1024
    bpm1 = 6

    bpm_data.read_tbt(home / "tbt_090513_215959.log", lin_opt)
    prt_FFT(cut, bpm_data.data[:, bpm1 - 1], window)

    window = 2
    cut = 0
    n_turn = 2048

    est_lin_opt.zero(len(bpm_data.loc))
    est_lin_opt.beta_pinger = np.array((6.92, 6.76), dtype=float)

    with open("tbt_optics.out", "w", encoding="utf-8") as outf:
        for file_name in (
            "tbt_090513_215619.log",
            "tbt_090513_215631.log",
            "tbt_090513_215652.log",
        ):
            bpm_data.read_tbt(home / file_name, lin_opt)
            get_nus(outf, cut, n_turn, window, bpm_data, lin_opt, est_lin_opt)

    est_lin_opt.get_stats(bpm_data, lin_opt)

    assert False

    cut = 10
    n_turn = 1024
    bpm1 = 5
    bpm2 = 6

    bpm_data.read_tbt(home / "tbt_090513_220010.log", lin_opt)
    ss_est(cut, n_turn, bpm1, bpm2, bpm_data, lin_opt)


def build_parser() -> argparse.ArgumentParser:
    parser = \
        argparse.ArgumentParser(
            description="Refactored turn-by-turn optics analysis")
    parser.add_argument(
        "home_dir", help="Directory containing linlat.out and TBT log files")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run(args.home_dir)


if __name__ == "__main__":
    main()
