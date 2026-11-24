#!/usr/bin/env python3
"""Extract day outlines from BbPN map with OpenCV; georef to obrys_pozaru.dat."""
from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(r"d:\wildfire")
KOD = ROOT / "kod"
MAP = str(ROOT / "img" / "biebrza_zasieg_pozaru_2020.jpg")
FINAL = np.loadtxt(KOD / "obrys_pozaru.dat")
UNIT_HA = 2.25
AX, BX, AY, BY = 111.554, 5.8286, 600.496, -5.4220


def poly_area(p):
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def close(p):
    return p if np.allclose(p[0], p[-1]) else np.vstack([p, p[0]])


def largest_contour(mask):
    cnts, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not cnts:
        return None
    c = max(cnts, key=cv2.contourArea)
    return c.reshape(-1, 2).astype(float)


def simplify(cnt, eps_frac=0.008):
    peri = cv2.arcLength(cnt.astype(np.float32).reshape(-1, 1, 2), True)
    approx = cv2.approxPolyDP(cnt.astype(np.float32).reshape(-1, 1, 2), eps_frac * peri, True)
    return approx.reshape(-1, 2)


def fit_affine(src, dst, n=64):
    def polar(p):
        p = p[:-1] if np.allclose(p[0], p[-1]) else p
        c = p.mean(0)
        a = np.arctan2(p[:, 1] - c[1], p[:, 0] - c[0])
        return p[np.argsort(a)]

    def resample(p, n):
        p = polar(p)
        nxt = np.vstack([p[1:], p[0]])
        seg = np.sqrt(((nxt - p) ** 2).sum(1))
        s = np.concatenate([[0], np.cumsum(seg)])
        t = np.linspace(0, s[-1], n, endpoint=False)
        x = np.interp(t, s, np.concatenate([p[:, 0], p[:1, 0]]))
        y = np.interp(t, s, np.concatenate([p[:, 1], p[:1, 1]]))
        return np.c_[x, y]

    s = resample(src, n)
    d = resample(dst, n)
    best, bestA = 1e99, None
    for sh in range(n):
        ss = np.roll(s, sh, 0)
        X = np.c_[ss, np.ones(n)]
        cx, *_ = np.linalg.lstsq(X, d[:, 0], rcond=None)
        cy, *_ = np.linalg.lstsq(X, d[:, 1], rcond=None)
        pred = np.c_[X @ cx, X @ cy]
        err = np.mean(((pred - d) ** 2).sum(1))
        if err < best:
            best, bestA = err, (cx, cy)
    print(f"affine RMSE={best**0.5:.2f}")
    return bestA


def apply(pts, A):
    X = np.c_[pts[:, 0], pts[:, 1], np.ones(len(pts))]
    return np.c_[X @ A[0], X @ A[1]]


def main():
    bgr = cv2.imread(MAP)
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB).astype(np.float32)
    h, w = rgb.shape[:2]
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    hsv = cv2.cvtColor(bgr, cv2.COLOR_BGR2HSV)

    # red-ish in HSV (wraps around 0)
    m1 = cv2.inRange(hsv, (0, 40, 40), (15, 255, 255))
    m2 = cv2.inRange(hsv, (160, 40, 40), (180, 255, 255))
    red = ((m1 > 0) | (m2 > 0)) & (R > G) & (R > B)
    # crop UI
    red[:, int(w * 0.78) :] = False
    red[: int(h * 0.09), :] = False
    red[int(h * 0.90) :, :] = False
    red[:, : int(w * 0.07)] = False
    red = cv2.morphologyEx(red.astype(np.uint8) * 255, cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8))
    red = cv2.morphologyEx(red, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))

    # brightness tiers inside red (V channel) — darker = earlier
    V = hsv[:, :, 2].astype(np.float32)
    vals = V[red > 0]
    print("V percentiles", np.percentile(vals, [15, 35, 55, 75, 90]))

    # cumulative masks by V threshold (darker first)
    def cum(pct):
        thr = np.percentile(vals, pct)
        m = (red > 0) & (V <= thr)
        m = cv2.morphologyEx(m.astype(np.uint8) * 255, cv2.MORPH_CLOSE, np.ones((7, 7), np.uint8))
        return m

    # tune pct so areas after georef are close — iterate
    m23 = red.copy()
    c23 = largest_contour(m23)
    c23s = simplify(c23, 0.004)
    A = fit_affine(c23s, FINAL)

    # search percentile cutoffs for target areas
    def area_for_mask(m):
        c = largest_contour(m)
        if c is None:
            return 0, None
        c = simplify(c, 0.006)
        dom = close(apply(c, A))
        return poly_area(dom) * UNIT_HA, dom

    # find pct for each target
    targets = [("20", 570), ("21", 1690), ("22", 4600)]
    pct_grid = np.linspace(8, 95, 88)
    chosen = {}
    for day, tgt in targets:
        best = None
        for pct in pct_grid:
            a, dom = area_for_mask(cum(pct))
            if dom is None:
                continue
            err = abs(a - tgt)
            if best is None or err < best[0]:
                best = (err, pct, a, dom)
        print(f"{day}: pct={best[1]:.1f} area={best[2]:.0f} (target {tgt})")
        chosen[day] = best[3]

    chosen["23"] = FINAL.copy()
    print(f"23: area={poly_area(FINAL)*UNIT_HA:.0f}")

    for day, poly in chosen.items():
        tgt = {"20": 570, "21": 1690, "22": 4600, "23": 5526}[day]
        a = poly_area(poly) * UNIT_HA
        # mild area scale about centroid if within reason
        if day != "23" and a > 50:
            s = float(np.clip((tgt / a) ** 0.5, 0.85, 1.2))
            c = poly[:-1].mean(0)
            poly = close((poly[:-1] - c) * s + c)
            a = poly_area(poly) * UNIT_HA
        np.savetxt(
            KOD / f"obrys_{day}.dat",
            poly,
            fmt="%.2f",
            header=f"Obrys {day}.04.2020 (~{a:.0f} ha, target {tgt}) z mapy BbPN (OpenCV); 1 j.=150 m",
        )
        print(f" wrote obrys_{day}.dat ({a:.0f} ha, n={len(poly)})")

    # preview
    fuel = Image.open(ROOT / "frames" / "fuel_map.png").convert("RGBA")
    ov = Image.new("RGBA", fuel.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(ov)

    def pix(p):
        return list(zip(AX + BX * p[:, 0], AY + BY * p[:, 1]))

    cols = {
        "22": (240, 140, 40, 50),
        "21": (220, 60, 0, 70),
        "20": (180, 0, 0, 110),
    }
    for day, fill in cols.items():
        p = np.loadtxt(KOD / f"obrys_{day}.dat")
        d.polygon(pix(p), fill=fill)
        d.line(pix(p), fill=(fill[0], fill[1], fill[2], 230), width=2)
    d.line(pix(FINAL), fill=(0, 0, 0, 220), width=2)
    Image.alpha_composite(fuel, ov).convert("RGB").save(ROOT / "frames" / "_debug_outlines_domain.png")
    print("preview ok")


if __name__ == "__main__":
    main()
