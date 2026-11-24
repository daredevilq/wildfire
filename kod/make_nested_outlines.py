#!/usr/bin/env python3
"""
Build nested historical outlines by clipping the hand-digitized 23.04 outline
along the NE→SW fire-growth axis, targeting reported areas (570/1690/4600/5526 ha).

Shapes follow the real final scar (obrys_pozaru.dat); cut position is fitted to area.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(r"d:\wildfire")
KOD = ROOT / "kod"
FINAL = np.loadtxt(KOD / "obrys_pozaru.dat")
UNIT_HA = 2.25

# Plot transform (gnuplot frames)
AX, BX, AY, BY = 111.554, 5.8286, 600.496, -5.4220


def poly_area(poly: np.ndarray) -> float:
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def close(poly: np.ndarray) -> np.ndarray:
    if not np.allclose(poly[0], poly[-1]):
        return np.vstack([poly, poly[0]])
    return poly


def clip_halfplane(poly: np.ndarray, origin: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Keep part of polygon where (p - origin)·normal >= 0 (Sutherland–Hodgman)."""
    poly = poly[:-1] if np.allclose(poly[0], poly[-1]) else poly
    nrm = normal / (np.linalg.norm(normal) + 1e-12)

    def inside(p):
        return np.dot(p - origin, nrm) >= -1e-9

    def intersect(p, q):
        dp = q - p
        denom = np.dot(dp, nrm)
        if abs(denom) < 1e-12:
            return q
        t = np.dot(origin - p, nrm) / denom
        return p + t * dp

    out = []
    for i in range(len(poly)):
        cur, prev = poly[i], poly[i - 1]
        cin, pin = inside(cur), inside(prev)
        if cin:
            if not pin:
                out.append(intersect(prev, cur))
            out.append(cur)
        elif pin:
            out.append(intersect(prev, cur))
    if len(out) < 3:
        return close(poly)
    return close(np.asarray(out, float))


def fit_cut(target_ha: float) -> np.ndarray:
    """
    Fire grows from NE tip toward SW. normal points toward NE (kept side).
    origin slides along axis until clipped area ≈ target_ha.
    """
    ne = np.array([71.4, 93.6])
    sw = np.array([12.9, 7.1])
    axis = ne - sw  # toward NE
    axis = axis / np.linalg.norm(axis)
    # normal for half-plane keep = toward NE = axis
    normal = axis

    # search origin = sw + t * (ne-sw), t in [0,1]
    best_poly, best_err = None, 1e99
    for t in np.linspace(0.05, 0.98, 200):
        origin = sw + t * (ne - sw)
        poly = clip_halfplane(FINAL, origin, normal)
        a = poly_area(poly) * UNIT_HA
        err = abs(a - target_ha)
        if err < best_err:
            best_err, best_poly = err, (poly, a, t)
    poly, a, t = best_poly
    print(f"  cut t={t:.3f} area={a:.0f} ha (target {target_ha}, err {abs(a-target_ha):.0f})")
    return poly


def to_pix(p):
    return list(zip(AX + BX * p[:, 0], AY + BY * p[:, 1]))


def main():
    specs = [
        ("20", 570),
        ("21", 1690),
        ("22", 4600),
        ("23", 5526),
    ]
    polys = {}
    for day, target in specs:
        print(day)
        if day == "23":
            poly = FINAL.copy()
            print(f"  final area={poly_area(poly)*UNIT_HA:.0f} ha")
        else:
            poly = fit_cut(target)
        # lightly densify for smoother dash
        polys[day] = poly
        hdr = (
            f"Obrys historyczny {day}.04.2020 (~{poly_area(poly)*UNIT_HA:.0f} ha, target {target}); "
            "wycinek obrysu 23.04 wzdluz osi NE-SW; 1 j.=150 m"
        )
        np.savetxt(KOD / f"obrys_{day}.dat", poly, fmt="%.2f", header=hdr)

    # preview on fuel_map
    fuel = Image.open(ROOT / "frames" / "fuel_map.png").convert("RGBA")
    ov = Image.new("RGBA", fuel.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(ov)
    style = {
        "20": ((180, 0, 0, 110), (180, 0, 0, 230)),
        "21": ((220, 60, 0, 70), (220, 60, 0, 220)),
        "22": ((240, 140, 40, 45), (240, 120, 0, 200)),
        "23": ((0, 0, 0, 0), (0, 0, 0, 220)),
    }
    for day in ("22", "21", "20", "23"):
        fill, line = style[day]
        pts = to_pix(polys[day])
        if fill[-1]:
            d.polygon(pts, fill=fill)
        d.line(pts, fill=line, width=2)
    Image.alpha_composite(fuel, ov).convert("RGB").save(ROOT / "frames" / "_debug_outlines_domain.png")
    print("saved frames/_debug_outlines_domain.png")


if __name__ == "__main__":
    main()
