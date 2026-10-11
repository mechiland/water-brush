"""Pigments and Kubelka-Munk optics.

A pigment carries absorption K and scattering S per RGB channel, per unit
layer thickness, plus the transfer parameters of Curtis et al. 1997 sec. 4.5:
settling speed (Curtis's density rho), staining power and granulation.

K and S are derived from two reflectances an artist can judge by eye, the
colour of a unit-thickness layer over white (Rw) and over black (Rb), using
the inversion in Curtis et al. 1997 sec. 5.1.
"""
from dataclasses import dataclass, field

import numpy as np


def ks_from_reflectance(Rw, Rb):
    """Curtis 1997 sec. 5.1: unit-thickness K, S from reflectance over white and black."""
    Rw = np.asarray(Rw, dtype=np.float64)
    Rb = np.asarray(Rb, dtype=np.float64)
    if np.any(Rb >= Rw) or np.any(Rb <= 0) or np.any(Rw >= 1):
        raise ValueError("need 0 < Rb < Rw < 1 in every channel")
    a = 0.5 * (Rw + (Rb - Rw + 1.0) / Rb)
    b = np.sqrt(a * a - 1.0)
    x = (b * b - (a - Rw) * (a - 1.0)) / (b * (1.0 - Rw))
    S = (1.0 / b) * 0.5 * np.log((x + 1.0) / (x - 1.0))  # arccoth
    K = S * (a - 1.0)
    return K, S


@dataclass(frozen=True)
class Pigment:
    name: str
    Rw: tuple            # reflectance of a unit layer over white
    Rb: tuple            # reflectance of a unit layer over black
    settle_mm_s: float = 3e-4  # Stokes settling speed; plays the role of Curtis rho
    staining: float = 2.0      # Curtis omega: resistance to lifting
    granulation: float = 0.1   # Curtis gamma: preference for paper valleys
    K: np.ndarray = field(init=False, repr=False, compare=False)
    S: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        K, S = ks_from_reflectance(self.Rw, self.Rb)
        object.__setattr__(self, "K", K)
        object.__setattr__(self, "S", S)


# Single pigment used for M1. Its colours are estimates of a burnt-sienna
# masstone; M3 replaces this with the full tube palette.
BURNT_SIENNA = Pigment(
    name="burnt_sienna",
    Rw=(0.60, 0.33, 0.20),
    Rb=(0.06, 0.04, 0.03),
    settle_mm_s=3e-4,
    staining=2.0,
    granulation=0.1,
)


def km_layer(K, S, x):
    """Reflectance R and transmittance T of a layer of thickness x (KM)."""
    S = np.maximum(S, 1e-6)
    a = 1.0 + K / S
    b = np.sqrt(a * a - 1.0)
    bSx = b * S * x
    sh = np.sinh(bSx)
    ch = np.cosh(bSx)
    c = a * sh + b * ch
    return sh / c, b / c


def km_composite(R1, T1, R2, T2):
    """Layer 1 on top of layer 2 (Kubelka 1954, Curtis eq. in sec. 5.2)."""
    denom = 1.0 - R1 * R2
    R = R1 + T1 * T1 * R2 / denom
    T = T1 * T2 / denom
    return R, T
