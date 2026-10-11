"""Render the M1 comparison image: renders/m1_flat_wash.png.

Left: the dried wash at 1000 px. Right top: the wet wash 60 s after laying
it, same crop. Right bottom: the dried top-left corner at 6x. Below: mean
deposit against distance from the edge.
"""
import pathlib
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from waterbrush.acceptance import edge_profile, edge_ratio, flat_wash  # noqa: E402
from waterbrush.render import render  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parents[1] / "renders" / "m1_flat_wash.png"


def profile_plot(dist, prof, ratio, size=(1000, 260)):
    W, H = size
    img = Image.new("RGB", size, (255, 255, 255))
    g = ImageDraw.Draw(img)
    x0, x1, y0, y1 = 70, W - 20, 20, H - 40
    top = prof.max() * 1.1
    g.line([(x0, y1), (x1, y1)], fill=(0, 0, 0))
    g.line([(x0, y0), (x0, y1)], fill=(0, 0, 0))
    xmax = dist.max()
    pts = [(x0 + (x1 - x0) * dd / xmax, y1 - (y1 - y0) * p / top) for dd, p in zip(dist, prof)]
    g.line(pts, fill=(160, 70, 30), width=2)
    for x, y in pts:
        g.ellipse([x - 2, y - 2, x + 2, y + 2], fill=(160, 70, 30))
    for mm in range(0, int(xmax) + 1, 2):
        x = x0 + (x1 - x0) * mm / xmax
        g.line([(x, y1), (x, y1 + 4)], fill=(0, 0, 0))
        g.text((x - 4, y1 + 6), str(mm), fill=(0, 0, 0))
    g.text((x0 + 10, y1 + 22), "distance from wash edge (mm)", fill=(0, 0, 0))
    g.text((4, y0), "deposit", fill=(0, 0, 0))
    g.text((x1 - 330, y0), f"edge (<=1 mm) / interior (>=10 mm) = {ratio:.2f}", fill=(0, 0, 0))
    return img


def main():
    sim, mask = flat_wash()
    sim.run(60.0)
    wet = render(sim, 1000)
    sim.run_until_dry()
    dry = render(sim, 1000)
    ratio = edge_ratio(sim.d[0], mask, sim.dx)
    dist, prof = edge_profile(sim.d[0], mask, sim.dx, max_mm=20.0)

    crop = (170, 70, 370, 200)     # top-left corner of the wash, px
    zoom = 3
    cw, ch = (crop[2] - crop[0]) * zoom, (crop[3] - crop[1]) * zoom
    wet_c = wet.crop(crop).resize((cw, ch), Image.NEAREST)
    dry_c = dry.crop(crop).resize((cw, ch), Image.NEAREST)

    plot = profile_plot(dist, prof, ratio, size=(dry.width, 260))
    out = Image.new("RGB", (dry.width + cw + 10, max(dry.height + 30 + plot.height, 2 * ch + 60)), (255, 255, 255))
    g = ImageDraw.Draw(out)
    out.paste(dry, (0, 20))
    g.text((4, 4), f"dried flat wash, cold press, t = {sim.t:.0f} s", fill=(0, 0, 0))
    out.paste(wet_c, (dry.width + 10, 20))
    g.text((dry.width + 14, 4), "wet, t = 60 s (corner, 3x)", fill=(0, 0, 0))
    out.paste(dry_c, (dry.width + 10, ch + 50))
    g.text((dry.width + 14, ch + 34), "dried (corner, 3x)", fill=(0, 0, 0))
    out.paste(plot, (0, dry.height + 30))
    OUT.parent.mkdir(exist_ok=True)
    out.save(OUT)
    print(f"wrote {OUT}  dry at t={sim.t:.0f}s  edge ratio {ratio:.3f}")


if __name__ == "__main__":
    main()
