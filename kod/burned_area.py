#!/usr/bin/env python3
# Powierzchnia wypalona w symulacji w kolejnych krokach (porownanie z pozarem Biebrzy 2020).
# Uruchamiac w katalogu z plikami fuel_*.data:  python3 burned_area.py > area.dat
# Punkt uznajemy za wypalony, jesli mial paliwo (> 0.2) i stracil co najmniej 40% paliwa.
import os

import numpy as np

H = 25.0  # godzin rzeczywistych na 1 s symulacji (hours_per_s w fire.hpp)
DT, LAST, EVERY = 1e-3, 6000, 100
UNIT_HA = 0.15 * 0.15 * 100  # 1 j.^2 = (150 m)^2 = 2.25 ha
OUTLINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "obrys_pozaru.dat")


def load(i):
    return np.loadtxt(f"fuel_{i}.data")


def smoothstep(a, b, t):
    h = np.clip((t - a) / (b - a), 0.0, 1.0)
    return h * h * (3 - 2 * h)


def forest(px, py):
    """Lasy jak w fire.hpp (funkcja forest): suma elips, 1 - las, 0 - brak lasu."""
    patches = [(36, 90, 18, 9), (12, 58, 8, 12)] + [(cx, cy, 5, 5) for cx, cy in
                                                    zip([20, 28, 40, 55, 62], [26, 38, 50, 65, 58])]
    f = np.zeros(px.shape)
    for x0, y0, a, b in patches:
        d = np.sqrt(((px - x0) / a) ** 2 + ((py - y0) / b) ** 2)
        f = np.maximum(f, 1 - smoothstep(0.8, 1.2, d))
    return f


def inside(px, py, poly):
    """Test punkt-w-wielokacie (ray casting) dla tablic punktow."""
    res = np.zeros(px.shape, dtype=bool)
    for (x1, y1), (x2, y2) in zip(poly[:-1], poly[1:]):
        cross = (y1 > py) != (y2 > py)
        xint = x1 + (py - y1) * (x2 - x1) / (y2 - y1 + 1e-300)
        res ^= cross & (px < xint)
    return res


d0 = load(0)
x, y, f0 = d0[:, 0], d0[:, 1], d0[:, 2]
n = int(round(np.sqrt(f0.size)))
cell = (100.0 / (n - 1)) ** 2 * UNIT_HA  # powierzchnia jednego punktu siatki wyjsciowej [ha]
is_forest = forest(x, y) > 0.5  # punkt lezy w placie lasu (polowa przejscia na brzegu elipsy)

print("# krok  t[s]  godz.  data  godzina  wypalone[ha]  w_lesie[ha]")
burned = None
for i in range(0, LAST + 1, EVERY):
    f = load(i)[:, 2]
    burned = (f0 > 0.2) & (f < 0.6 * f0)
    h = i * DT * H
    total = 20 + h  # t = 0 <-> 19.04 20:00
    m = int(round(60 * total))  # minuty od 19.04 00:00 (round - bez bledow zaokraglen)
    day, hour, minute = 19 + m // 1440, (m // 60) % 24, m % 60
    print(f"{i:5d} {i * DT:5.1f} {h:6.1f}  {day:02d}.04 {hour:02d}:{minute:02d}  "
          f"{burned.sum() * cell:7.0f}  {(burned & is_forest).sum() * cell:7.0f}")

# zgodnosc ksztaltu z rzeczywistym obrysem (23.04.2020) dla ostatniego kroku
poly = np.loadtxt(OUTLINE)
real = inside(x, y, poly)
common = (burned & real).sum()
print(f"# obrys rzeczywisty: {real.sum() * cell:.0f} ha, czesc wspolna: {common * cell:.0f} ha, "
      f"IoU = {common / (burned | real).sum():.2f}, "
      f"pokrycie obrysu = {common / real.sum():.2f}, trafnosc symulacji = {common / burned.sum():.2f}")
