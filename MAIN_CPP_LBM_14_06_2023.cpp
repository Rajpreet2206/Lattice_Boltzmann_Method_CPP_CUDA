
#include <iostream>
#include <cmath>
#include<chrono>

int main() {
    int nx = 501, ny = 201;
    int i, j, a, a1, ts, ia, ja;
    float** f[9], ** ft[9]; // store current and a temporary distribution
    int isn[201][21];
    float wt[9], ux, uy, rho, term1, term2, u2;
    int ex[9], ey[9], kb[9];
    float source[9];
    int time = 2000;
    float dpdx = 1.0e-5;
    float tau = 0.8;
    float rho0 = 1.0;
    float visc, H;
    float feq[9];
    float rho_av;
    float ux_exact[21];
    int kpor;
    float delta = 0.5;
    float y, y2;
    float c1, c2;
    float l2;

    visc = (tau - 0.5) / 3.0;
    H = ny - 1 - 2.0 * delta;

    ex[0] = 0;
    ey[0] = 0;
    ex[1] = 1;
    ey[1] = 0;
    ex[2] = 0;
    ey[2] = 1;
    ex[3] = -1;
    ey[3] = 0;
    ex[4] = 0;
    ey[4] = -1;
    ex[5] = 1;
    ey[5] = 1;
    ex[6] = -1;
    ey[6] = 1;
    ex[7] = -1;
    ey[7] = -1;
    ex[8] = 1;
    ey[8] = -1;

    for (a = 0; a < 9; a++) {
        f[a] = new float*[nx];
        ft[a] = new float*[nx];
        for (i = 0; i < nx; i++) {
            f[a][i] = new float[ny];
            ft[a][i] = new float[ny];
        }
    }

    for (a = 0; a < 9; a++) {
        if (a == 0) { wt[a] = 4.0 / 9.0; }
        if (a >= 1 && a <= 4) { wt[a] = 1.0 / 9.0; }
        if (a >= 5 && a <= 8) { wt[a] = 1.0 / 36.0; }
    }

    for (a = 0; a < 9; a++) {
        for (a1 = a; a1 < 9; a1++) {
            if ((ex[a] + ex[a1] == 0 && ey[a] + ey[a1] == 0)) {
                kb[a] = a1;
                kb[a1] = a;
            }
        }
    }
    for (a = 0; a < 9; a++) {
        std::cout << "a = " << a << ", kb = " << kb[a] << std::endl;
    }

    // Initialization
    auto Tstart = std::chrono::high_resolution_clock::now();
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            isn[i][j] = 0;
            if (j == 0 || j == ny - 1) {
                isn[i][j] = 1;
            }
            for (a = 0; a < 9; a++) {
                f[a][i][j] = wt[a] * rho0; // this corresponds to zero velocity
            }
        }
    }
    auto Tend = std::chrono::high_resolution_clock::now();
    auto time_duration = std::chrono::duration_cast<std::chrono::milliseconds>(Tend - Tstart).count();
    std::cout<< time_duration << " milliseconds " << std::endl;

auto TstartSL = std::chrono::high_resolution_clock::now();
    for (ts = 1; ts <= time; ts++) {
        rho_av = 0.0;
        kpor = 0;
        for (i = 0; i < nx; i++) {
            for (j = 1; j < ny - 1; j++) {
                ux = 0.0;
                uy = 0.0;
                rho = 0.0;
                for (a = 0; a < 9; a++) {
                    rho += f[a][i][j];
                    ux += f[a][i][j];
                    uy += f[a][i][j];
                }
                rho_av += rho;
                kpor += 1;
                ux += dpdx / 2.0;
                ux /= rho;
                uy /= rho;
                u2 = ux * ux + uy * uy;
                // ux += tau * dpdx
                for (a = 0; a < 9; a++) {
                    term1 = ux * ex[a] + uy * ey[a];
                    term2 = term1 * term1;
                    source[a] = (1.0 - 0.5 / tau) * wt[a] * (3.0 * (ex[a] - ux) + 9.0 * (ex[a] * ux + ey[a] * uy) * ex[a]) * dpdx;
                    feq[a] = wt[a] * rho * (1.0 + 3.0 * term1 + 4.5 * term2 - 1.5 * u2);
                    ft[a][i][j] = f[a][i][j] - (f[a][i][j] - feq[a]) / tau + source[a]; // collision step
                }
            }
        }
/*
        if (ts % 100 == 0) {
            std::cout << "ts = " << ts << "\t rho_av = " << rho_av / kpor << std::endl;
        }
*/

        // Streaming of the particle dist. after collision
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny - 1; j++) {
                for (a = 0; a < 9; a++) {
                    ia = i + ex[a];
                    ja = j + ey[a];
                    if (ia < 0) { ia = nx - 1; }
                    if (ia > nx - 1) { ia = 0; }
                    f[a][ia][ja] = ft[a][i][j];
                }
            }
        }
        // Boundary condition, especially bounceback
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                if (isn[i][j] == 0) {
                    for (a = 0; a < 9; a++) {
                        ia = i - ex[a];
                        ja = j - ey[a];
                        if (ia < 0) { ia = nx - 1; }
                        if (ia > nx - 1) { ia = 0; }
                        if (isn[ia][ja] == 1) {
                            f[a][i][j] = f[kb[a]][ia][ja];
                        }
                    }
                }
            }
        }
    }
auto TendSL = std::chrono::high_resolution_clock::now();
auto time_durationSL = std::chrono::duration_cast<std::chrono::milliseconds>(TendSL - TstartSL).count();
std::cout<< time_durationSL << " milliseconds " << std::endl;

    for (j = 0; j < ny; j++) {
        y = (j - delta) / H;
        y2 = y * y;
        ux_exact[j] = 0.5 * dpdx * H * H * (y - y2) / visc;
        if (j == 0 || j == ny - 1) {
            ux_exact[j] = 0;
        }
    }

    i = nx - 1;
    l2 = 0.0;
    c1 = 0.0;
    c2 = 0.0;
    for (j = 0; j < ny; j++) {
        rho = 0.0;
        ux = 0.0;
        uy = 0.0;
        if (isn[i][j] == 0) {
            for (a = 0; a < 9; a++) {
                rho += f[a][i][j];
                ux += f[a][i][j] * ex[a];
                uy += f[a][i][j] * ey[a];
            }
            ux /= rho;
            uy /= rho;
        }
    }
    c1 += ux_exact[j] * ux_exact[j];
    c2 += (ux_exact[j] - ux) * (ux_exact[j] - ux);
    std::cout << nx - 1 << " " << j << " " << ux << " " << ux_exact[j] << " " << uy << " " << rho << std::endl;
    std::cout << "c1 = " << c1 << " c2 = " << c2 << std::endl;
    l2 = std::pow((c2 / c1), 0.5);
    std::cout << "l2 = " << l2 << std::endl;

    // Clean up memory
    for (a = 0; a < 9; a++) {
        for (i = 0; i < nx; i++) {
            delete[] f[a][i];
            delete[] ft[a][i];
        }
        delete[] f[a];
        delete[] ft[a];
    }

    return 0;
}