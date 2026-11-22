# Wykresy i animacje do sprawozdania. Uruchamiac w katalogu z plikami out_*.data, fuel_*.data, wind_*.data:
#   python3 /code/examples/fire/burned_area.py > area.dat
#   gnuplot /code/examples/fire/plot.gp
# Wyniki: frames/*.png, anim.gif (temperatura), fuel_anim.gif (paliwo), wind_anim.gif (wiatr)

H = 25.0          # godzin rzeczywistych na 1 s symulacji (hours_per_s w fire.hpp)
dt = 1e-3         # krok czasowy
last = 6000       # ostatni krok
every_n = 100     # co ile krokow sa zapisane pliki
outline = "/code/examples/fire/obrys_pozaru.dat"   # rzeczywisty obrys pozaru (23.04.2020)

# data i godzina rzeczywista dla kroku i (t = 0 <-> 19.04.2020 20:00)
hrs(i) = 20 + i*dt*H
mins(i) = int(60 * hrs(i) + 0.5)   # minuty od 19.04 00:00 (zaokraglenie - bez bledow typu 02:29)
date(i) = sprintf("%02d.04 %02d:%02d", 19 + mins(i) / 1440, (mins(i) / 60) % 24, mins(i) % 60)
label(i) = sprintf("krok %d (t = %.1f s, ok. %s)", i, i*dt, date(i))

set size ratio -1
set xrange [0:100]
set yrange [0:100]
set xlabel "x [j.]  (1 j. = 150 m, wschod)"
set ylabel "y [j.]  (polnoc)"
set tics out
set key bottom right opaque

temp_palette = "set palette defined (300 '#2c7bb6', 550 '#abd9e9', 800 '#ffffbf', 1100 '#fdae61', 1500 '#d7191c'); set cbrange [300:1500]; set cblabel 'Temperatura [K]'"
fuel_palette = "set palette defined (0 '#3b2a20', 0.2 '#6e5040', 0.4 '#a08c3c', 0.6 '#c5d63a', 1 '#1a7a2a'); set cbrange [0:1]; set cblabel 'Paliwo (ulamek)'"

system("mkdir -p frames")
set terminal pngcairo size 860,700 font ",11"

# --- mapa paliwa i stan poczatkowy ---
@fuel_palette
set output "frames/fuel_map.png"
set title "Mapa paliwa fuel\\_init(x, y)"
plot "fuel_0.data" using 1:2:3 with image notitle, \
     outline with lines lw 2 dt 2 lc rgb "#d01010" title "obrys rzeczywisty (23.04)"

@temp_palette
set output "frames/init_state.png"
set title "Stan poczatkowy init\\_state(x, y)"
plot "out_0.data" using 1:2:3 with image notitle

# --- klatki temperatury, paliwa i wiatru ---
do for [i in "0 500 1000 1500 2000 2500 3000 3400 4000 6000"] {
    n = int(i)
    @temp_palette
    set output sprintf("frames/temp_%05d.png", n)
    set title "Temperatura - ".label(n)
    plot sprintf("out_%d.data", n) using 1:2:3 with image notitle

    @fuel_palette
    set output sprintf("frames/fuel_%05d.png", n)
    set title "Paliwo - ".label(n)
    plot sprintf("fuel_%d.data", n) using 1:2:3 with image notitle, \
         outline with lines lw 2 dt 2 lc rgb "#d01010" title "obrys rzeczywisty (23.04)"

    @fuel_palette
    set output sprintf("frames/wind_%05d.png", n)
    set title "Wiatr b = (bx, by) - ".label(n)
    plot "fuel_0.data" using 1:2:3 with image notitle, \
         sprintf("wind_%d.data", n) using 1:2:($3/10):($4/10) with vectors head filled size 1,20 lc rgb "#1040c0" notitle
}

# --- powierzchnia wypalona: symulacja vs. pozar Biebrzy (mapa BbPN) ---
$real << EOD
13.5  570
37.3 1690
61.2 4600
85.6 5526
EOD
if (system("test -f area.dat && echo 1") eq "1") {
    reset
    set terminal pngcairo size 800,500 font ",11"
    set output "frames/area.png"
    set title "Powierzchnia wypalona - symulacja i pozar w BbPN (2020)"
    set xlabel "czas od zaplonu (19.04 ok. 20:00) [h]"
    set ylabel "powierzchnia [ha]"
    set grid
    set key bottom right
    plot "area.dat" using 3:6 with lines lw 2 lc rgb "#d01010" title "symulacja", \
         $real using 1:2 with linespoints pt 7 ps 1.5 lw 1 dt 2 lc rgb "#202020" title "pozar BbPN 2020 (mapa zasiegu)"
    set size ratio -1
    set xrange [0:100]
    set yrange [0:100]
    unset grid
}

# --- animacje ---
# "set terminal ... animate" trzeba powtorzyc przed kazdym plikiem: zmiana "set output" konczy animacje
# i kolejne klatki bylyby zapisane jako osobne, nieanimowane obrazki
anim_term = "set terminal gif animate delay 10 loop 0 size 860,700 font ',11'"
set xlabel "x [j.]  (1 j. = 150 m, wschod)"
set ylabel "y [j.]  (polnoc)"

@anim_term
@temp_palette
set output "anim.gif"
do for [n = 0:last:every_n] {
    set title "Temperatura - ".label(n)
    plot sprintf("out_%d.data", n) using 1:2:3 with image notitle
}

@anim_term
@fuel_palette
set output "fuel_anim.gif"
do for [n = 0:last:every_n] {
    set title "Paliwo - ".label(n)
    plot sprintf("fuel_%d.data", n) using 1:2:3 with image notitle, \
         outline with lines lw 2 dt 2 lc rgb "#d01010" title "obrys rzeczywisty (23.04)"
}

@anim_term
set output "wind_anim.gif"
do for [n = 0:last:every_n] {
    set title "Wiatr b = (bx, by) - ".label(n)
    plot "fuel_0.data" using 1:2:3 with image notitle, \
         sprintf("wind_%d.data", n) using 1:2:($3/10):($4/10) with vectors head filled size 1,20 lc rgb "#1040c0" notitle
}
unset output
