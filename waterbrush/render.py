"""Render the deposition layer (and any pigment still in the water) with KM optics."""
import numpy as np
from PIL import Image

from .pigment import km_layer


PAPER_WHITE = np.array([0.93, 0.92, 0.88])   # reflectance of bare cotton paper (estimate)


def _resample(a, H, W):
    """Bilinear resample of a 2-D field to H x W (cell centres aligned)."""
    h0, w0 = a.shape
    y = (np.arange(H) + 0.5) * h0 / H - 0.5
    x = (np.arange(W) + 0.5) * w0 / W - 0.5
    y0 = np.clip(np.floor(y).astype(int), 0, h0 - 1)
    x0 = np.clip(np.floor(x).astype(int), 0, w0 - 1)
    y1 = np.clip(y0 + 1, 0, h0 - 1)
    x1 = np.clip(x0 + 1, 0, w0 - 1)
    fy = np.clip(y - y0, 0, 1)[:, None]
    fx = np.clip(x - x0, 0, 1)[None, :]
    a00 = a[y0][:, x0]
    a01 = a[y0][:, x1]
    a10 = a[y1][:, x0]
    a11 = a[y1][:, x1]
    return (a00 * (1 - fx) + a01 * fx) * (1 - fy) + (a10 * (1 - fx) + a11 * fx) * fy


def reflectance(sim, width_px=1000, include_wet=True):
    """Linear RGB reflectance image, shape (H_px, W_px, 3).

    KM is evaluated per simulation cell, then the reflectance image is
    resampled to the output size, so a pixel straddling an edge gets the
    average colour of the cells it covers.
    """
    shade = 1.0 - 0.06 * sim.paper.relief_mm / 0.08 * (1.0 - sim.h)
    R = PAPER_WHITE[None, None, :] * shade[:, :, None]
    # Pigment layers, bottom first: deposited pigment, then suspended pigment.
    layers = [sim.d] + ([sim.G] if include_wet else [])
    for amounts in layers:
        for k, pg in enumerate(sim.pigments):
            x = amounts[k]
            if not np.any(x > 0):
                continue
            Rl, Tl = km_layer(pg.K[None, None, :], pg.S[None, None, :], x[:, :, None])
            R = Rl + Tl * Tl * R / (1.0 - Rl * R)   # this layer over what is below
    W = int(width_px)
    H = int(round(sim.H * W / sim.W))
    return np.stack([_resample(R[:, :, ch], H, W) for ch in range(3)], axis=-1)


def to_srgb8(R):
    R = np.clip(R, 0.0, 1.0)
    srgb = np.where(R <= 0.0031308, 12.92 * R, 1.055 * np.power(R, 1 / 2.4) - 0.055)
    return (np.clip(srgb, 0, 1) * 255 + 0.5).astype(np.uint8)


def render(sim, width_px=1000, include_wet=True):
    return Image.fromarray(to_srgb8(reflectance(sim, width_px, include_wet)))
