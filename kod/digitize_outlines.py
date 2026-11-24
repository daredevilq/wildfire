#!/usr/bin/env python3
"""
Digitize cumulative historical burn outlines from BbPN extent map,
georeferenced to domain using obrys_pozaru.dat (23.04).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(r"d:\wildfire")
KOD = ROOT / "kod"
MAP = ROOT / "img" / "biebrza_zasieg_pozaru_2020.jpg"
FINAL = np.loadtxt(KOD / "obrys_pozaru.dat")
UNIT_HA = 2.25


def poly_area(poly: np.ndarray) -> float:
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def order_contour(ys: np.ndarray, xs: np.ndarray) -> np.ndarray:
    """Order boundary pixels into a closed polyline (nearest-neighbor walk)."""
    pts = np.c_[xs.astype(float), ys.astype(float)]
    # start at topmost
    i0 = int(np.argmin(pts[:, 1] - 0.01 * pts[:, 0]))
    used = np.zeros(len(pts), dtype=bool)
    order = [i0]
    used[i0] = True
    for _ in range(len(pts) - 1):
        last = pts[order[-1]]
        d = np.sum((pts - last) ** 2, axis=1)
        d[used] = np.inf
        j = int(np.argmin(d))
        if not np.isfinite(d[j]):
            break
        # jump break -> stop this component
        if d[j] > 40**2:
            break
        order.append(j)
        used[j] = True
    ring = pts[order]
    if len(ring) > 2:
        ring = np.vstack([ring, ring[0]])
    return ring


def boundary_pixels(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    from numpy.lib.stride_tricks import sliding_window_view

    m = mask.astype(np.uint8)
    # erode manually
    pad = np.pad(m, 1, constant_values=0)
    win = sliding_window_view(pad, (3, 3))
    interior = win.min(axis=(-1, -2))
    edge = m.astype(bool) & ~interior.astype(bool)
    ys, xs = np.where(edge)
    return ys, xs


def largest_component(mask: np.ndarray) -> np.ndarray:
    """Simple flood-fill largest True component."""
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    best = None
    best_n = 0
    ys, xs = np.where(mask)
    for y0, x0 in zip(ys, xs):
        if seen[y0, x0]:
            continue
        stack = [(y0, x0)]
        seen[y0, x0] = True
        cells = []
        while stack:
            y, x = stack.pop()
            cells.append((y, x))
            for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                yy, xx = y + dy, x + dx
                if 0 <= yy < h and 0 <= xx < w and mask[yy, xx] and not seen[yy, xx]:
                    seen[yy, xx] = True
                    stack.append((yy, xx))
        if len(cells) > best_n:
            best_n = len(cells)
            best = cells
    out = np.zeros_like(mask)
    if best:
        for y, x in best:
            out[y, x] = True
    return out


def subsample_ring(ring: np.ndarray, n: int = 40) -> np.ndarray:
    """Return closed ring with exactly n+1 points (last=first), or open n points if n requested for open."""
    ring = np.asarray(ring, float)
    if len(ring) < 3:
        raise ValueError("ring too short")
    if np.allclose(ring[0], ring[-1]):
        open_r = ring[:-1]
    else:
        open_r = ring
    # arc-length resample of open part, then close
    nxt = np.vstack([open_r[1:], open_r[0]])
    seg = np.sqrt(((nxt - open_r) ** 2).sum(axis=1))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(s[-1])
    if total < 1e-6:
        raise ValueError("zero-length ring")
    # n open samples around the loop
    targets = np.linspace(0, total, n, endpoint=False)
    # s has len n_open+1; interpolate with wrap by extending
    xs = np.interp(targets, s, np.concatenate([open_r[:, 0], open_r[:1, 0]]))
    ys = np.interp(targets, s, np.concatenate([open_r[:, 1], open_r[:1, 1]]))
    out = np.c_[xs, ys]
    return np.vstack([out, out[0]])


def fit_affine(src: np.ndarray, dst: np.ndarray):
    """src (map px x,y), dst (domain) -> affine domain = f(px).

    Uses polar-sorted, equal-length resampling + circular shift search.
    """
    def polar_sort(p):
        p = np.asarray(p, float)
        if np.allclose(p[0], p[-1]):
            p = p[:-1]
        c = p.mean(axis=0)
        ang = np.arctan2(p[:, 1] - c[1], p[:, 0] - c[0])
        return p[np.argsort(ang)]

    n = 48
    s = subsample_ring(polar_sort(src), n)[:-1]
    d = subsample_ring(polar_sort(dst), n)[:-1]
    assert s.shape == (n, 2) and d.shape == (n, 2)

    best_err, best_A = np.inf, None
    for shift in range(n):
        ss = np.roll(s, shift, axis=0)
        X = np.c_[ss[:, 0], ss[:, 1], np.ones(n)]
        cx, *_ = np.linalg.lstsq(X, d[:, 0], rcond=None)
        cy, *_ = np.linalg.lstsq(X, d[:, 1], rcond=None)
        pred = np.c_[X @ cx, X @ cy]
        err = float(np.mean(np.sum((pred - d) ** 2, axis=1)))
        if err < best_err:
            best_err = err
            best_A = (cx, cy)
    print(f"affine fit RMSE={best_err**0.5:.3f} domain units")
    return best_A


def apply_affine(pts_xy: np.ndarray, A) -> np.ndarray:
    cx, cy = A
    X = np.c_[pts_xy[:, 0], pts_xy[:, 1], np.ones(len(pts_xy))]
    return np.c_[X @ cx, X @ cy]


def mask_to_domain_poly(mask: np.ndarray, A, max_pts: int = 36) -> np.ndarray:
    mask = largest_component(mask)
    ys, xs = boundary_pixels(mask)
    if len(xs) < 20:
        # fallback: all True pixels bbox corners — should not happen
        ys, xs = np.where(mask)
    # subsample boundary before ordering for speed
    step = max(len(xs) // 800, 1)
    ring_px = order_contour(ys[::step], xs[::step])
    ring_px = subsample_ring(ring_px, max_pts)
    dom = apply_affine(ring_px, A)
    # close
    if not np.allclose(dom[0], dom[-1]):
        dom = np.vstack([dom, dom[0]])
    return dom


def main():
    im = np.asarray(Image.open(MAP).convert("RGB")).astype(np.float32)
    h, w = im.shape[:2]
    R, G, B = im[:, :, 0], im[:, :, 1], im[:, :, 2]

    # burn scar: reddish, exclude UI (legend right, title top, scale bottom)
    reddish = (R > 120) & (R > G + 25) & (R > B + 25) & (G < 180)
    burn = reddish.copy()
    burn[:, int(w * 0.78) :] = False
    burn[: int(h * 0.10), :] = False
    burn[int(h * 0.90) :, :] = False
    burn[:, : int(w * 0.08)] = False
    burn = largest_component(burn)
    print("burn pixels", int(burn.sum()))

    # Separate by brightness matching legend: darker = earlier
    # Sample legend colors roughly from right side legend swatches if possible,
    # else use R+G value tiers inside burn.
    val = (0.6 * R + 0.3 * G + 0.1 * B)
    v = val[burn]
    # cumulative masks: day20 = darkest quantile, day21 = darker half, etc.
    # Use fixed thresholds tuned to map shades (dark red / mid / light / pale)
    # Inspect percentiles
    print("val percentiles", np.percentile(v, [10, 25, 50, 75, 90]))

    # thresholds from visual legend: dark maroon, red, light red, pale
    # cumulative:
    m20 = burn & (val < np.percentile(v, 22))  # earliest core
    m21 = burn & (val < np.percentile(v, 48))
    m22 = burn & (val < np.percentile(v, 78))
    m23 = burn

    # Grow earlier days to be nested (morphological): ensure m20 ⊂ m21 ⊂ m22 ⊂ m23
    # Already true by construction if percentiles nest.

    # Fit affine using day23 outer contour vs FINAL outline
    ys, xs = boundary_pixels(m23)
    step = max(len(xs) // 1000, 1)
    ring23 = subsample_ring(order_contour(ys[::step], xs[::step]), 48)
    A = fit_affine(ring23, FINAL)

    # Debug overlay on map
    dbg = Image.fromarray(im.astype(np.uint8))
    draw = ImageDraw.Draw(dbg, "RGBA")
    for mask, col in [(m20, (255, 0, 0, 80)), (m21, (255, 128, 0, 50)), (m22, (255, 255, 0, 40))]:
        ys, xs = np.where(mask)
        for y, x in zip(ys[::20], xs[::20]):
            draw.point((x, y), fill=col)
    dbg.save(ROOT / "frames" / "_debug_burn_tiers.png")

    targets = {"20": (m20, 570), "21": (m21, 1690), "22": (m22, 4600), "23": (m23, 5526)}
    for day, (mask, target_ha) in targets.items():
        # For day 23 use the hand-digitized final outline (already calibrated)
        if day == "23":
            poly = FINAL.copy()
        else:
            poly = mask_to_domain_poly(mask, A, max_pts=32)
            # scale about centroid to better match reported area (keep shape)
            a = poly_area(poly) * UNIT_HA
            if a > 10:
                s = (target_ha / a) ** 0.5
                # limit correction — don't destroy georef
                s = float(np.clip(s, 0.75, 1.35))
                c = poly[:-1].mean(axis=0)
                poly = np.vstack([(poly[:-1] - c) * s + c, (poly[0] - c) * s + c])
                # re-close
                poly[-1] = poly[0]
        a = poly_area(poly) * UNIT_HA
        path = KOD / f"obrys_{day}.dat"
        hdr = (
            f"Obrys historyczny {day}.04.2020 (~{a:.0f} ha, target {target_ha}); "
            "z mapy BbPN, georef do domeny; 1 j.=150 m"
        )
        np.savetxt(path, poly, fmt="%.2f", header=hdr)
        print(f"{day}.04: {a:.0f} ha (target {target_ha}), n={len(poly)} -> {path.name}")

    # Also write preview in domain on fuel_map
    try:
        fuel = Image.open(ROOT / "frames" / "fuel_map.png").convert("RGBA")
        # use same plot transform as rebuild_cmp2
        ax, bx, ay, by = 111.554, 5.8286, 600.496, -5.4220

        def to_pix(p):
            return list(zip(ax + bx * p[:, 0], ay + by * p[:, 1]))

        ov = Image.new("RGBA", fuel.size, (0, 0, 0, 0))
        d = ImageDraw.Draw(ov)
        colors = {
            "20": (180, 0, 0, 100),
            "21": (220, 80, 0, 70),
            "22": (240, 160, 40, 50),
            "23": (0, 0, 0, 0),
        }
        for day, col in colors.items():
            p = np.loadtxt(KOD / f"obrys_{day}.dat")
            pts = to_pix(p)
            if day != "23":
                d.polygon(pts, fill=col)
            d.line(pts, fill=(col[0], col[1], col[2], 220), width=2)
        # final outline black dashed-ish
        d.line(to_pix(FINAL), fill=(0, 0, 0, 200), width=2)
        Image.alpha_composite(fuel, ov).convert("RGB").save(ROOT / "frames" / "_debug_outlines_domain.png")
        print("saved frames/_debug_outlines_domain.png")
    except Exception as e:
        print("preview failed", e)


if __name__ == "__main__":
    main()
