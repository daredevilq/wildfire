// SPDX-FileCopyrightText: 2015 - 2024 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef FIRE_FIRE_HPP
#define FIRE_FIRE_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

#include "ads/executor/galois.hpp"
#include "ads/output_manager.hpp"
#include "ads/simulation.hpp"
#include "ads/util.hpp"

inline double falloff(double r, double R, double t) {
    if (t < r)
        return 1.0;
    if (t > R)
        return 0.0;
    double h = (t - r) / (R - r);
    return std::pow((h - 1) * (h + 1), 2);
}

// Zarzewie o srodku (x0, y0)
inline double bump(double r, double R, double x, double y, double x0 = 50, double y0 = 50) {
    double dx = x - x0;
    double dy = y - y0;
    double t = std::sqrt(dx * dx + dy * dy) / 100;
    return falloff(r / 200, R / 200, t);
}

// Gladkie przejscie 0 -> 1 dla t w [a, b] (dziala tez dla a > b)
inline double smoothstep(double a, double b, double t) {
    double h = std::clamp((t - a) / (b - a), 0.0, 1.0);
    return h * h * (3 - 2 * h);
}

inline double mix(double a, double b, double h) { return (1 - h) * a + h * b; }

namespace ads::problems {

class fire : public simulation_2d {
private:
    using Base = simulation_2d;
    vector_type u, u_prev;

    vector_type fuel, fuel_prev;

    // double C1 = 1;
    // double C2 = 1;
    // double C3 = 1;
    // double C4 = 1;
    // double C5 = 1;
    // double C6 = 1;
    // double T0 = 1;
    // double bx = 40;
    // double by = 20;
    // Wiatr: funkcje bx(x, y, t), by(x, y, t) - patrz nizej

    // Skale scenariusza "Pozar Biebrzanskiego Parku Narodowego, 19-26.04.2020":
    // 1 jednostka domeny = 150 m (domena 100 x 100 = 15 km x 15 km, os y na polnoc),
    // 1 s symulacji = hours_per_s godzin rzeczywistych, t = 0 <-> 19.04.2020 ok. godz. 20:00
    double hours_per_s = 25.0;

    double ch = 1.0;
    double Ar = 5.7e-5;
    double rho = 1.293;
    double cp = 1.0;
    double cw = 0.5;
    double kappa = 0.3;
    double xi = 2e-2;
    double sigma = 5.67e-8;
    double eps = 0.05; // "percentage error"
    double delta_x = 3.5e-2 / eps;
    double delta_z = 1.5 * eps;
    double T0 = 300;
    double Tcomb = 1200;
    double Tig = 800;
    double Ta = 300; // ???

    double M = 2;
    double M1 = 1;

    galois_executor executor;
    output_manager<2> output;

public:
    explicit fire(const config_2d& config, int threads)
    : Base{config}
    , u{shape()}
    , u_prev{shape()}
    , fuel{shape()}
    , fuel_prev{shape()}
    , executor{threads}
    , output{x.B, y.B, 300} { }

    // Stan poczatkowy: dwa male zarzewia (wypalanie traw) w polnocnej czesci Basenu Srodkowego Biebrzy
    double init_state(double x, double y) {
        double r = 2;   // pelna temperatura w promieniu 1 j. (150 m)
        double R = 16;  // brak podgrzania dalej niz 8 j. (1200 m), lagodne zbocze - stabilnosc
        double main_spot = bump(r, R, x, y, 72, 82);
        double second_spot = bump(r, 0.75 * R, x, y, 67, 73);
        return T0 + Tcomb * std::max(main_spot, second_spot);
    }

    // Odleglosc ze znakiem od koryta Biebrzy (prosta (13, 7) -> (82, 65)); > 0 na SE (pola uprawne)
    static double river_dist(double x, double y) { return 0.644 * (x - 13) - 0.766 * (y - 7); }

    // Odleglosc ze znakiem od skraju suchej czesci doliny (prosta (12, 24) -> (50, 82)); > 0 na NW
    static double edge_dist(double x, double y) { return -0.836 * (x - 12) + 0.548 * (y - 24); }

    // Plat lasu w ksztalcie elipsy o srodku (x0, y0) i polosiach a, b: 1 - las, 0 - brak lasu
    static double ellipse(double x, double y, double x0, double y0, double a, double b) {
        double d = std::sqrt(std::pow((x - x0) / a, 2) + std::pow((y - y0) / b, 2));
        return 1 - smoothstep(0.8, 1.2, d);
    }

    // Lasy - oddzielne platy (las nie jest ciagly)
    static double forest(double x, double y) {
        double f = std::max(ellipse(x, y, 36, 90, 18, 9),   // Czerwone Bagno i Grzedy
                            ellipse(x, y, 12, 58, 8, 12));  // Brzeziny Kapickie
        const double px[] = {20, 28, 40, 55, 62};           // male platy lasu (brzezina, ols) w dolinie
        const double py[] = {26, 38, 50, 65, 58};
        for (int i = 0; i < 5; ++i) {
            f = std::max(f, ellipse(x, y, px[i], py[i], 5, 5));
        }
        return f;
    }

    // Mapa paliwa (zakres 0..1, 1 - latwopalne, 0 - nie ma co sie palic)
    double fuel_map(double x, double y) const {
        double f = 1.0;  // suche trzcinowiska, turzycowiska i laki (susza 2020)
        f = mix(f, 0.15, smoothstep(-2, 2, edge_dist(x, y)));    // NW: wilgotne laki i bagna - nie pala sie
        f = mix(f, 0.15, smoothstep(92, 96, y));                 // N: mokre olsy
        f = mix(f, 0.4, forest(x, y));                           // lasy - malo suchej biomasy w runie
        f = mix(f, 0.1, smoothstep(82, 86, x));                  // pola i enklawy na E
        f = mix(f, 0.1, smoothstep(22, 18, x + y));              // pola przy Goniadzu (SW)
        f = mix(f, 0.1, smoothstep(-1, 1, river_dist(x, y)));    // za rzeka: pola, wsie
        f *= 1 - std::exp(-std::pow(river_dist(x, y) / 1.2, 2));  // koryto Biebrzy - brak paliwa
        double ex = x - 76, ey = y - 76;
        f = mix(f, 0.1, 1 - smoothstep(3, 5, std::sqrt(ex * ex + ey * ey)));  // enklawa rolna
        return f;
    }

    // Kierunek wiatru (kat wektora predkosci w stopniach, 0 = na wschod, 90 = na polnoc)
    // i jego predkosc w zaleznosci od czasu rzeczywistego h [godz. od 19.04 20:00]
    double wind_angle(double h) const {
        double a = -80;                                // wiatr z NNW -> ogien idzie na S wzdluz rzeki (20.04)
        a = mix(a, -135, smoothstep(12, 24, h));       // wiatr z NE  -> ogien idzie na SW (21-22.04)
        a = mix(a, -180, smoothstep(58, 64, h));       // 22.04 rano: wiatr z E -> ogien na W
        a = mix(a, -270, smoothstep(64, 74, h));       // 22-23.04: wiatr z S -> ogien skreca na N
        return a;
    }

    double wind_speed(double h) const {
        double w = 25;
        w = mix(w, 40, smoothstep(20, 30, h));   // silny wiatr 21-22.04
        w = mix(w, 20, smoothstep(56, 66, h));   // oslabienie 22.04
        w = mix(w, 5, smoothstep(80, 90, h));    // 23.04 - prawie cisza, pozar opanowany
        return w;
    }

    // Wiatr bx(x, y, t), by(x, y, t); nad lasem wiatr jest slabszy
    double wind_factor(double x, double y) const { return 1 - 0.4 * forest(x, y); }

    double bx(double x, double y, double t) const {
        double h = t * hours_per_s;
        return wind_factor(x, y) * wind_speed(h) * std::cos(wind_angle(h) * M_PI / 180);
    }

    double by(double x, double y, double t) const {
        double h = t * hours_per_s;
        return wind_factor(x, y) * wind_speed(h) * std::sin(wind_angle(h) * M_PI / 180);
    }

private:
    void before() override {
        prepare_matrices();

        auto init = [this](double x, double y) { return init_state(x, y); };
        projection(u, init);
        solve(u);
        output.to_file(u, "out_0.data");

        // Mapa paliwa
        auto fuel_init = [this](double x, double y) { return fuel_map(x, y); };
        projection(fuel, fuel_init);
        solve(fuel);
        output.to_file(fuel, "fuel_0.data");
        save_wind(0, 0.0);
    }

    void before_step(int /*iter*/, double /*t*/) override {
        using std::swap;
        swap(u, u_prev);
        swap(fuel, fuel_prev);
    }

    void step(int /*iter*/, double t) override {
        compute_rhs(t);
        solve(u);
        solve(fuel);
        apply_bc(u);
    }

    // Warunek brzegowy T = T0 na brzegu domeny (daleko od pozaru - temperatura otoczenia).
    // Bez niego wiatr "wpadajacy" do domeny przez brzeg powoduje niestabilnosc schematu jawnego.
    void apply_bc(vector_type& T) const {
        int nx = x.dofs();
        int ny = y.dofs();
        for (int i = 0; i < nx; ++i) {
            T(i, 0) = T0;
            T(i, ny - 1) = T0;
        }
        for (int j = 0; j < ny; ++j) {
            T(0, j) = T0;
            T(nx - 1, j) = T0;
        }
    }

    void compute_rhs(double t) {
        auto& rhs = u;
        auto& rhs_fuel = fuel;

        zero(rhs);
        zero(rhs_fuel);
        executor.for_each(elements(), [&](index_type e) {
            auto U = element_rhs();
            auto F = element_rhs();

            double J = jacobian(e);
            for (auto q : quad_points()) {
                double w = weight(q);
                // Wiatr zalezny od polozenia i czasu
                auto p = point(e, q);
                double x = p[0];
                double y = p[1];
                double wx = bx(x, y, t);
                double wy = by(x, y, t);

                value_type u = eval_fun(u_prev, e, q);
                value_type fuel = eval_fun(fuel_prev, e, q);
                // Zabezpieczenie numeryczne - przy silnym wietrze oscylacje schematu jawnego
                // moga dac T < 0 K, a wtedy czlony ~T^3 i ~T^4 "uciekaja" do -inf (NaN).
                // W czlonach nieliniowych uzywamy wiec T = max(u, 1 K).
                double T = std::max(u.val, 1.0);

                for (auto a : dofs_on_element(e)) {
                    auto aa = dof_global_to_local(e, a);
                    value_type v = eval_basis(e, q, a);
                    double inv = 1.0 / (rho * cp);

                    double delta = u.val > Tig && fuel.val > 0.2 ? 1.0 : 0.0;
                    double hc = -70; // enthalpy
                    double r = delta * Ar * T * std::exp(-Ta / T);
                    double Rc = - 1e4 * rho * ch * hc * M / M1 * r;
                    double Qw = - rho * cw * (wx * u.dx + wy * u.dy);
                    double qc = - kappa * grad_dot(u, v);
                    double qd = 0.0; // omitted
                    double qr = -4 * sigma * eps * delta_x * std::pow(T, 3) * grad_dot(u, v);
                    double Qconv = xi * (T0 - u.val);
                    double Qrz = sigma * eps / delta_z * (std::pow(T0, 4) - std::pow(T, 4));

                    double val = (Rc + Qw + Qconv + Qrz) * v.val + qc + qd + qr;
                    U(aa[0], aa[1]) += (u.val * v.val + steps.dt * inv * val) * w * J;

                    double fval = - delta * 3e2 * r * fuel.val * v.val;
                    F(aa[0], aa[1]) += (fuel.val * v.val + steps.dt * fval) * w * J;
                }
            }
            executor.synchronized([&] { update_global_rhs(rhs, U, e); });
            executor.synchronized([&] { update_global_rhs(rhs_fuel, F, e); });
        });
    }

    void after_step(int iter, double t) override {
        auto i = iter  +1;
        // Zapis wynikow co 100 krokow + zapis pola wiatru
        if (i % 100 == 0) {
            std::cout << "Step " << i << std::endl;
            output.to_file(u, "out_%d.data", i);
            output.to_file(fuel, "fuel_%d.data", i);
            save_wind(i, t + steps.dt);
        }
    }

    // Zapis pola wiatru na siatce 21 x 21 punktow: "x y bx by" (do gnuplota: with vectors)
    void save_wind(int i, double t) const {
        char name[64];
        std::snprintf(name, sizeof(name), "wind_%d.data", i);
        std::ofstream os{name};
        for (int ix = 0; ix <= 20; ++ix) {
            for (int iy = 0; iy <= 20; ++iy) {
                double x = 5.0 * ix;
                double y = 5.0 * iy;
                os << x << ' ' << y << ' ' << bx(x, y, t) << ' ' << by(x, y, t) << '\n';
            }
            os << '\n';
        }
    }

    // double forcing(point_type x, double /*t*/) const {
    //     return 0;
    // }
};

}  // namespace ads::problems

#endif  // FIRE_FIRE_HPP
