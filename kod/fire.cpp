// SPDX-FileCopyrightText: 2015 - 2024 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "fire.hpp"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: fire <N> <p> <threads>" << std::endl;
        std::exit(1);
    }
    int n = std::atoi(argv[1]);
    int p = std::atoi(argv[2]);
    int threads = std::atoi(argv[3]);

    // p=2, n=200
    ads::dim_config dim{p, n, 0.0, 100.0};
    // 6000 krokow: 6 s = 150 h przy hours_per_s = 25, tj. do 26.04 (koniec akcji gasniczej)
    ads::timesteps_config steps{6000, 1e-3};
    int ders = 1;

    ads::config_2d c{dim, dim, steps, ders};
    ads::problems::fire sim{c, threads};
    sim.run();
}
