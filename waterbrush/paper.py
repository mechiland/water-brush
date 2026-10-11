"""Paper: surface height field and absorbency parameters.

The height field h in [0, 1] drives three things in the simulation, all taken
from Curtis et al. 1997 (sec. 4.2-4.5): water collects in valleys, granulating
pigment settles in valleys, and capillary capacity follows h. It is generated
procedurally from a seed so replay is exact.
"""
from dataclasses import dataclass

import numpy as np

# Texture presets. tooth_mm is the feature size of the paper grain, relief_mm
# the peak-to-valley depth of the surface. See notes/research.md (paper).
TEXTURES = {
    "hot": {"tooth_mm": 0.6, "relief_mm": 0.01},
    "cold": {"tooth_mm": 1.5, "relief_mm": 0.05},
    "rough": {"tooth_mm": 2.5, "relief_mm": 0.10},
}


def _smooth_noise(rng, shape, sigma_cells):
    """White noise low-passed with a Gaussian of the given width (in cells)."""
    noise = rng.standard_normal(shape)
    if sigma_cells <= 0.3:
        return noise
    H, W = shape
    fy = np.fft.fftfreq(H)[:, None]
    fx = np.fft.rfftfreq(W)[None, :]
    kernel = np.exp(-2.0 * (np.pi * sigma_cells) ** 2 * (fx * fx + fy * fy))
    return np.fft.irfft2(np.fft.rfft2(noise) * kernel, s=shape)


@dataclass(frozen=True)
class Paper:
    texture: str = "cold"   # "hot", "cold" or "rough"
    sizing: float = 0.8     # 0 = waterleaf (unsized), 1 = hard sized
    weight_gsm: float = 300.0

    def __post_init__(self):
        if self.texture not in TEXTURES:
            raise ValueError(f"unknown paper texture {self.texture!r}")
        if not 0.0 <= self.sizing <= 1.0:
            raise ValueError("sizing must be in [0, 1]")
        if self.weight_gsm <= 0:
            raise ValueError("weight_gsm must be positive")

    @property
    def relief_mm(self):
        return TEXTURES[self.texture]["relief_mm"]

    @property
    def thickness_mm(self):
        # Cotton watercolour paper: 300 gsm is about 0.5 mm thick.
        return 0.5 * self.weight_gsm / 300.0

    @property
    def capacity_mm(self):
        """Water the paper can hold per unit area, as a depth in mm.

        Thickness times an effective porosity of 0.3 (estimate).
        """
        return 0.3 * self.thickness_mm

    @property
    def sorptivity(self):
        """Lucas-Washburn sorptivity A, in mm/s^0.5: absorbed depth = A*sqrt(t).

        Calibrated to the Cobb60 test (water taken up in 60 s; 1 g/m^2 is
        1 um of water depth): about 15 g/m^2 for hard-sized paper and 150
        g/m^2 for unsized paper (estimates, see notes/research.md).
        """
        cobb60_mm = 0.001 * (15.0 + 135.0 * (1.0 - self.sizing) ** 2)
        return cobb60_mm / np.sqrt(60.0)

    def height_field(self, shape, dx_mm, seed):
        """Surface height h in [0, 1] for a grid of the given shape."""
        rng = np.random.default_rng([seed, 0x9A9E2])
        tooth = TEXTURES[self.texture]["tooth_mm"]
        # Two octaves: the tooth itself plus a finer fibre structure.
        h = _smooth_noise(rng, shape, 0.5 * tooth / dx_mm)
        h = h / (h.std() + 1e-12)
        fine = _smooth_noise(rng, shape, 0.25 * tooth / dx_mm)
        h = h + 0.35 * fine / (fine.std() + 1e-12)
        lo, hi = np.percentile(h, [0.5, 99.5])
        return np.clip((h - lo) / (hi - lo), 0.0, 1.0)
