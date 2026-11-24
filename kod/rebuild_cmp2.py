#!/usr/bin/env python3
"""Rebuild comparison frames: remove old red by color, add day-specific historical area."""
from __future__ import annotations

import subprocess
from io import BytesIO
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(r"d:\wildfire")
FRAMES = ROOT / "frames"
KOD = ROOT / "kod"
H, DT = 25.0, 1e-3
# refined plot mapping from red-cloud fit
A_X, B_X = 111.554, 5.8286
A_Y, B_Y = 600.496, -5.4220


def hours(step):
    return step * DT * H


def date_str(step):
    total = 20 + hours(step)
    m = int(round(60 * total))
    return f"{19 + m // 1440:02d}.04 {(m // 60) % 24:02d}:{m % 60:02d}"


def label(step):
    return f"krok {step} (t = {hours(step):.1f} h, {date_str(step)})"


def font(size=15):
    for n in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\calibri.ttf"):
        if Path(n).exists():
            return ImageFont.truetype(n, size)
    return ImageFont.load_default()


def load_orig(step):
    data = subprocess.check_output(["git", "show", f"HEAD:frames/fuel_{step:05d}.png"], cwd=str(ROOT))
    return Image.open(BytesIO(data)).convert("RGB")


def to_pix(poly):
    return list(zip((A_X + B_X * poly[:, 0]).tolist(), (A_Y + B_Y * poly[:, 1]).tolist()))


def remove_red_and_key(im: Image.Image) -> Image.Image:
    arr = np.asarray(im).astype(np.float64)
    h, w = arr.shape[:2]
    R, G, B = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    # strict + soft red (dashed anti-alias)
    core = (R > 160) & (G < 90) & (B < 90) & (R > G + 70) & (R > B + 70)
    soft = (R > 140) & (R > G + 35) & (R > B + 35) & (G < 160) & (B < 140) & (R > 1.15 * G)
    mask = core | soft
    # keep colorbar (right) and title
    mask[:, 720:] = False
    mask[:50, :] = False
    # dilate
    from scipy import ndimage  # may not exist

    try:
        mask = ndimage.binary_dilation(mask, iterations=2)
    except Exception:
        m2 = mask.copy()
        ys, xs = np.where(mask)
        for dy in range(-2, 3):
            for dx in range(-2, 3):
                m2[np.clip(ys + dy, 0, h - 1), np.clip(xs + dx, 0, w - 1)] = True
        mask = m2

    for _ in range(6):
        pad = np.pad(arr, ((1, 1), (1, 1), (0, 0)), mode="edge")
        acc = np.zeros_like(arr)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dx == dy == 0:
                    continue
                acc += pad[1 + dy : h + 1 + dy, 1 + dx : w + 1 + dx]
        arr[mask] = (acc / 8.0)[mask]

    # always paint over gnuplot key area (bottom-right of plot) with local terrain color
    sample = arr[500:520, 480:680]
    fill = np.median(sample.reshape(-1, 3), axis=0)
    arr[530:625, 460:710] = fill
    return Image.fromarray(np.clip(arr, 0, 255).astype(np.uint8))


def cover_title(im, title):
    out = im.copy()
    d = ImageDraw.Draw(out)
    w, _ = out.size
    d.rectangle([0, 0, w, 48], fill=(255, 255, 255))
    f = font(15)
    bb = d.textbbox((0, 0), title, font=f)
    d.text(((w - (bb[2] - bb[0])) / 2, 14), title, fill=(0, 0, 0), font=f)
    return out


def add_historical(im, poly, legend):
    out = im.convert("RGBA")
    overlay = Image.new("RGBA", out.size, (0, 0, 0, 0))
    od = ImageDraw.Draw(overlay)
    pts = to_pix(poly)
    od.polygon(pts, fill=(220, 30, 30, 55), outline=(200, 0, 0, 200))
    # dashed outline on top
    base = Image.alpha_composite(out, overlay).convert("RGB")
    d = ImageDraw.Draw(base)
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        segs = max(int(np.hypot(x1 - x0, y1 - y0) / 6), 1)
        for s in range(segs):
            if s % 2:
                continue
            t0, t1 = s / segs, min((s + 1) / segs, 1)
            d.line(
                [(x0 + (x1 - x0) * t0, y0 + (y1 - y0) * t0), (x0 + (x1 - x0) * t1, y0 + (y1 - y0) * t1)],
                fill=(200, 0, 0),
                width=3,
            )
    f = font(13)
    bx0, by0 = 455, 555
    bb = d.textbbox((0, 0), legend, font=f)
    tw, th = bb[2] - bb[0], bb[3] - bb[1]
    d.rectangle([bx0 - 8, by0 - 6, bx0 + tw + 40, by0 + th + 10], fill=(255, 255, 255), outline=(0, 0, 0))
    d.rectangle([bx0, by0 + 2, bx0 + 18, by0 + th - 2], fill=(220, 30, 30, 255) if False else (220, 80, 80))
    d.text((bx0 + 26, by0), legend, fill=(0, 0, 0), font=f)
    return base


def main():
    specs = [
        (500, "20", "historycznie spalone do 20.04"),
        (1500, "21", "historycznie spalone do 21.04"),
        (2500, "22", "historycznie spalone do 22.04"),
        (3400, "23", "historycznie spalone do 23.04"),
    ]
    for step, day, leg in specs:
        im = load_orig(step)
        im = remove_red_and_key(im)
        im = cover_title(im, "Paliwo - " + label(step))
        poly = np.loadtxt(KOD / f"obrys_{day}.dat")
        im = add_historical(im, poly, leg)
        im.save(FRAMES / f"fuel_{step:05d}.png")
        im.save(FRAMES / f"fuel_cmp_{step:05d}.png")
        # verify no leftover large red outside day poly bbox
        arr = np.asarray(im)
        R, G, B = arr[:, :, 0].astype(int), arr[:, :, 1].astype(int), arr[:, :, 2].astype(int)
        red = (R > 160) & (G < 100) & (B < 100) & (R > G + 60)
        red[:, 720:] = False
        red[:50, :] = False
        ys, xs = np.where(red)
        print(f"{step}: red_px={len(xs)} y[{ys.min() if len(ys) else 0}:{ys.max() if len(ys) else 0}]")


if __name__ == "__main__":
    main()
