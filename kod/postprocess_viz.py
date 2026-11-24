#!/usr/bin/env python3
"""
Postprocess frames/GIFs for the Biebrza fire report:
1) titles: t in hours (not simulation seconds)
2) temp + terrain + wind composite animation
3) fuel comparison frames with day-specific historical burn outlines
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(r"d:\wildfire")
FRAMES = ROOT / "frames"
KOD = ROOT / "kod"

H = 25.0
DT = 1e-3

# domain -> pixel (fitted to gnuplot 860x700 frames with colorbar)
A_X, B_X = 111.554, 5.8286
A_Y, B_Y = 600.496, -5.4220

# plot box (approx) — used to limit overlays
PLOT = (110, 55, 695, 605)  # left, top, right, bottom


def hours(step: int) -> float:
    return step * DT * H


def date_str(step: int) -> str:
    total = 20 + hours(step)  # t=0 <-> 19.04 20:00
    m = int(round(60 * total))
    day, hour, minute = 19 + m // 1440, (m // 60) % 24, m % 60
    return f"{day:02d}.04 {hour:02d}:{minute:02d}"


def label(step: int) -> str:
    return f"krok {step} (t = {hours(step):.1f} h, {date_str(step)})"


def font(size: int = 16):
    for name in (
        r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibri.ttf",
        r"C:\Windows\Fonts\DejaVuSans.ttf",
    ):
        p = Path(name)
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def to_pix(poly: np.ndarray) -> list[tuple[float, float]]:
    xs = A_X + B_X * poly[:, 0]
    ys = A_Y + B_Y * poly[:, 1]
    return list(zip(xs.tolist(), ys.tolist()))


def cover_title(im: Image.Image, new_title: str) -> Image.Image:
    """Cover top title band and redraw title (gnuplot-style)."""
    out = im.copy()
    draw = ImageDraw.Draw(out)
    w, _ = out.size
    # title sits in white margin above plot
    draw.rectangle([0, 0, w, 48], fill=(255, 255, 255))
    f = font(15)
    # center-ish title
    bbox = draw.textbbox((0, 0), new_title, font=f)
    tw = bbox[2] - bbox[0]
    draw.text(((w - tw) / 2, 14), new_title, fill=(0, 0, 0), font=f)
    return out


def parse_step_from_name(name: str) -> int | None:
    m = re.search(r"_(\d{5})\.png$", name)
    return int(m.group(1)) if m else None


def fix_static_frames():
    mapping = {
        "temp_": "Temperatura - ",
        "fuel_": "Paliwo - ",
        "wind_": "Wiatr b = (bx, by) - ",
    }
    for path in sorted(FRAMES.glob("*.png")):
        prefix = None
        for p, title in mapping.items():
            if path.name.startswith(p):
                prefix = title
                break
        if prefix is None:
            continue
        step = parse_step_from_name(path.name)
        if step is None:
            continue
        im = Image.open(path).convert("RGB")
        im = cover_title(im, prefix + label(step))
        im.save(path)
        print("title", path.name)


def fix_gif(path: Path, title_prefix: str, every: int = 100):
    im = Image.open(path)
    frames = []
    n = getattr(im, "n_frames", 1)
    for i in range(n):
        im.seek(i)
        step = i * every
        fr = cover_title(im.convert("RGB"), title_prefix + label(step))
        frames.append(fr)
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        duration=im.info.get("duration", 100),
        loop=0,
    )
    print(f"gif titles {path.name}: {n} frames")


def load_outline(day: str) -> np.ndarray:
    return np.loadtxt(KOD / f"obrys_{day}.dat")


def erase_old_outline(im: Image.Image) -> Image.Image:
    """Paint over the known 23.04 outline path (more reliable than color thresholding)."""
    arr = np.asarray(im).copy().astype(np.float64)
    h, w = arr.shape[:2]
    poly = np.loadtxt(KOD / "obrys_pozaru.dat")
    pts = to_pix(poly)

    mask = np.zeros((h, w), dtype=bool)
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        n = max(int(np.hypot(x1 - x0, y1 - y0) * 2), 1)
        xs = np.round(np.linspace(x0, x1, n)).astype(int)
        ys = np.round(np.linspace(y0, y1, n)).astype(int)
        for x, y in zip(xs, ys):
            y0_, y1_ = max(0, y - 4), min(h, y + 5)
            x0_, x1_ = max(0, x - 4), min(w, x + 5)
            mask[y0_:y1_, x0_:x1_] = True

    # blur-fill: replace masked pixels by local mean of non-masked neighbors
    for _ in range(4):
        pad = np.pad(arr, ((1, 1), (1, 1), (0, 0)), mode="edge")
        acc = np.zeros_like(arr)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dx == 0 and dy == 0:
                    continue
                acc += pad[1 + dy : h + 1 + dy, 1 + dx : w + 1 + dx]
        avg = acc / 8.0
        arr[mask] = avg[mask]
    return Image.fromarray(np.clip(arr, 0, 255).astype(np.uint8))


def draw_outline(im: Image.Image, poly: np.ndarray, legend: str) -> Image.Image:
    out = im.copy()
    draw = ImageDraw.Draw(out)
    pts = to_pix(poly)
    # dashed polyline
    for i in range(len(pts) - 1):
        x0, y0 = pts[i]
        x1, y1 = pts[i + 1]
        # dash
        segs = max(int(np.hypot(x1 - x0, y1 - y0) / 8), 1)
        for s in range(segs):
            if s % 2:
                continue
            t0, t1 = s / segs, min((s + 1) / segs, 1.0)
            draw.line(
                [
                    (x0 + (x1 - x0) * t0, y0 + (y1 - y0) * t0),
                    (x0 + (x1 - x0) * t1, y0 + (y1 - y0) * t1),
                ],
                fill=(208, 16, 16),
                width=2,
            )
    # legend box bottom-right inside plot
    f = font(12)
    text = legend
    bx0, by0 = 430, 560
    bbox = draw.textbbox((0, 0), text, font=f)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    pad = 6
    draw.rectangle(
        [bx0 - pad, by0 - pad, bx0 + tw + 28, by0 + th + pad],
        fill=(255, 255, 255),
        outline=(0, 0, 0),
    )
    draw.line([bx0, by0 + th // 2, bx0 + 18, by0 + th // 2], fill=(208, 16, 16), width=2)
    draw.text((bx0 + 24, by0 - 1), text, fill=(0, 0, 0), font=f)
    return out


def load_original_fuel_frame(step: int) -> Image.Image:
    """Load pristine fuel frame from git HEAD (before postprocess overlays)."""
    import subprocess

    data = subprocess.check_output(
        ["git", "show", f"HEAD:frames/fuel_{step:05d}.png"],
        cwd=str(ROOT),
    )
    from io import BytesIO

    return Image.open(BytesIO(data)).convert("RGB")


def make_comparison_frames():
    """4 fuel frames with historical burned area for that day."""
    specs = [
        (500, "20", "obrys rzeczywisty (20.04)"),
        (1500, "21", "obrys rzeczywisty (21.04)"),
        (2500, "22", "obrys rzeczywisty (22.04)"),
        (3400, "23", "obrys rzeczywisty (23.04)"),
    ]
    for step, day, leg in specs:
        src = FRAMES / f"fuel_{step:05d}.png"
        im = load_original_fuel_frame(step)
        im = erase_old_outline(im)
        # blank old gnuplot key (bottom-right of plot) before drawing ours
        draw = ImageDraw.Draw(im)
        sample = np.asarray(im)[505:520, 450:650]
        fill = tuple(int(x) for x in np.median(sample.reshape(-1, 3), axis=0))
        draw.rectangle([450, 545, 700, 620], fill=fill)
        im = cover_title(im, "Paliwo - " + label(step))
        im = draw_outline(im, load_outline(day), leg)
        out = FRAMES / f"fuel_cmp_{step:05d}.png"
        im.save(src)
        im.save(out)
        print("comparison", out.name)


def extract_wind_arrows(wind: Image.Image, fuel: Image.Image) -> np.ndarray:
    """Mask of wind-arrow pixels: blue-ish difference vs evolving fuel frame."""
    w = np.asarray(wind.convert("RGB"), dtype=np.int16)
    f = np.asarray(fuel.convert("RGB"), dtype=np.int16)
    diff = np.linalg.norm(w.astype(np.float32) - f.astype(np.float32), axis=2)
    B = w[:, :, 2]
    R = w[:, :, 0]
    G = w[:, :, 1]
    # blue arrows on wind frames
    arrows = (B > R + 20) & (B > G + 10) & (B > 80) & (diff > 25)
    arrows[:, 720:] = False
    arrows[:50, :] = False
    return arrows


def _strip_fuel_overlay(fuel_rgb: np.ndarray) -> np.ndarray:
    """Remove dashed red outline + gnuplot key from a fuel frame (keep fuel colors)."""
    arr = fuel_rgb.astype(np.float64)
    h, w = arr.shape[:2]
    R, G, B = arr[:, :, 0], arr[:, :, 1], arr[:, :, 2]
    mask = (R > 140) & (R > G + 35) & (R > B + 35) & (G < 160)
    mask[:, 720:] = False
    mask[:50, :] = False
    # dilate
    m2 = mask.copy()
    ys, xs = np.where(mask)
    for dy in range(-2, 3):
        for dx in range(-2, 3):
            m2[np.clip(ys + dy, 0, h - 1), np.clip(xs + dx, 0, w - 1)] = True
    mask = m2
    for _ in range(4):
        pad = np.pad(arr, ((1, 1), (1, 1), (0, 0)), mode="edge")
        acc = np.zeros_like(arr)
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dx == dy == 0:
                    continue
                acc += pad[1 + dy : h + 1 + dy, 1 + dx : w + 1 + dx]
        arr[mask] = (acc / 8.0)[mask]
    # cover key box
    sample = arr[500:520, 480:680]
    fill = np.median(sample.reshape(-1, 3), axis=0)
    arr[530:625, 460:710] = fill
    return arr


def composite_temp_fuel_wind_frame(
    temp: Image.Image, fuel: Image.Image, wind: Image.Image
) -> Image.Image:
    """Evolving fuel (depleting) + wind arrows + hot temperature overlay."""
    t = np.asarray(temp.convert("RGB"), dtype=np.float32)
    out = _strip_fuel_overlay(np.asarray(fuel.convert("RGB")))

    L, T, R, Btm = PLOT
    warr = np.asarray(wind.convert("RGB"))
    Ww = warr.astype(np.int16)
    blue = (Ww[:, :, 2] > Ww[:, :, 0] + 25) & (Ww[:, :, 2] > Ww[:, :, 1] + 15) & (Ww[:, :, 2] > 90)
    blue[:, 720:] = False
    blue[:50, :] = False
    out[blue] = warr[blue]

    region_t = t[T:Btm, L:R]
    region_f = out[T:Btm, L:R]
    ambient = np.array([44.0, 123.0, 182.0])
    dist = np.linalg.norm(region_t - ambient, axis=2)
    hot = dist > 55
    alpha = np.clip((dist - 55) / 80.0, 0.0, 0.85) * hot
    a = alpha[..., None]
    out[T:Btm, L:R] = (1 - a) * region_f + a * region_t

    out[:, 720:] = t[:, 720:]
    return Image.fromarray(np.clip(out, 0, 255).astype(np.uint8))


def make_temp_terrain_wind_anim():
    """Composite anim.gif (T) + fuel_anim.gif (depleting fuel) + wind arrows."""
    # Prefer pristine GIFs from git if current ones were title-patched
    import subprocess
    from io import BytesIO

    def open_gif(name: str) -> Image.Image:
        path = ROOT / name
        try:
            data = subprocess.check_output(["git", "show", f"HEAD:{name}"], cwd=str(ROOT))
            return Image.open(BytesIO(data))
        except Exception:
            return Image.open(path)

    temp_gif = open_gif("anim.gif")
    fuel_gif = open_gif("fuel_anim.gif")
    wind_gif = open_gif("wind_anim.gif")
    n = min(
        getattr(temp_gif, "n_frames", 1),
        getattr(fuel_gif, "n_frames", 1),
        getattr(wind_gif, "n_frames", 1),
    )
    frames = []
    for i in range(n):
        temp_gif.seek(i)
        fuel_gif.seek(i)
        wind_gif.seek(i)
        step = i * 100
        fr = composite_temp_fuel_wind_frame(
            temp_gif.convert("RGB"), fuel_gif.convert("RGB"), wind_gif.convert("RGB")
        )
        # strip fuel outline/key noise lightly? keep as-is; title overwrite
        fr = cover_title(fr, "Temperatura + paliwo + wiatr - " + label(step))
        frames.append(fr)
        if step in {0, 500, 1000, 1500, 2000, 2500, 3000, 3400, 4000, 6000}:
            fr.save(FRAMES / f"temp_terrain_wind_{step:05d}.png")
        if i % 15 == 0:
            print(f"  composite frame {i}/{n}")
    out = ROOT / "temp_terrain_wind_anim.gif"
    frames[0].save(
        out,
        save_all=True,
        append_images=frames[1:],
        duration=100,
        loop=0,
    )
    print(f"wrote {out.name} ({n} frames)")


def main():
    fix_static_frames()
    fix_gif(ROOT / "anim.gif", "Temperatura - ")
    fix_gif(ROOT / "fuel_anim.gif", "Paliwo - ")
    fix_gif(ROOT / "wind_anim.gif", "Wiatr b = (bx, by) - ")
    make_temp_terrain_wind_anim()
    # porównanie obszarów (historyczne obrysy dzienne) — osobny skrypt
    import runpy

    runpy.run_path(str(KOD / "rebuild_cmp2.py"), run_name="__main__")


if __name__ == "__main__":
    main()
