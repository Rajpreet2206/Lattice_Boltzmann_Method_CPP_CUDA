#include <iostream>
#include <cmath>

const int N = 50;            // Grid size
const double L = 1.0;        // Length of the channel
const double U0 = 0.1;       // Maximum velocity
const double nu = 0.1;       // Viscosity
const double rho0 = 1.0;     // Initial density or the reference density

const double dt = 1.0;       // Time step
const double dx = L / N;     // Spatial step

const int maxIter = 100000;    // Maximum number of iterations
const double convergence = 1e-5; // Convergence criteria

// Lattice velocity vectors
const int Q = 9;
const int cx[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
const int cy[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
const double weights[Q] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };//Center Velocity[0,], Axis-Aligned Velocities [1,2,3,4]

double tau = 3.0 * nu + 0.5;
double omega = 1.0 / tau;
//double u_conv = U0 * dx / maxVel;
//dt = dx * dx / (u_conv * u_conv);

// 2D arrays to store the distributions
double f[N][N][Q];   // Distribution functions at t
double f_next[N][N][Q];  // Distribution functions at t+1

// Function to initialize the distributions
void initialize() {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double u = U0 * 6.0 * (j * (N - j)) / (N * N);
            double rho = rho0;

            for (int k = 0; k < Q; ++k) {
                double cu = 3.0 * (cx[k] * u);
                double uu = 3.0 * (cx[k] * u) * (cx[k] * u);
                f[i][j][k] = weights[k] * rho * (1.0 + cu + 0.5 * uu - 1.5 * u * u);//Equillibria Discrete Velocities
            }
        }
    }
}

void collide_stream() {
    // Collision step
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double rho = 0.0, ux = 0.0, uy = 0.0;

            for (int k = 0; k < Q; ++k) {
                rho += f[i][j][k];
                ux += cx[k] * f[i][j][k];
                uy += cy[k] * f[i][j][k];
            }

            ux /= rho;
            uy /= rho;

            for (int k = 0; k < Q; ++k) {
                double cu = 3.0 * (cx[k] * ux + cy[k] * uy);
                double uu = 3.0 * (cx[k] * ux + cy[k] * uy) * (cx[k] * ux + cy[k] * uy);
                f_next[i][j][k] = f[i][j][k] - (f[i][j][k] - weights[k] * rho * (1.0 + cu + 0.5 * uu - 1.5 * (ux * ux + uy * uy))) / tau;
            }
        }
    }

    // Stream step
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < Q; ++k) {
                int inew = (i + cx[k] + N) % N;
                int jnew = (j + cy[k] + N) % N;
                f[inew][jnew][k] = f_next[i][j][k];
            }
        }
    }
}

// Function to calculate the maximum velocity in the channel
double maxVelocity() {
    double maxVel = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double ux = 0.0, uy = 0.0;

            for (int k = 0; k < Q; ++k) {
                ux += cx[k] * f[i][j][k];
                uy += cy[k] * f[i][j][k];
            }

            double vel = std::sqrt(ux * ux + uy * uy);
            if (vel > maxVel) {
                maxVel = vel;
            }
        } 
    }
    return maxVel;
}


int main(){
    initialize();

    int iter = 0;
    double residual = convergence + 1.0;
    while (iter < maxIter && residual > convergence) {
        double maxVel = maxVelocity();
        double tau = 3.0 * nu + 0.5;
        double omega = 1.0 / tau;
        double u_conv = U0 * dx / maxVel;
        const int dt = dx * dx / (u_conv * u_conv);

        collide_stream();

        double sum = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < Q; ++k) {
                    sum += std::abs(f[i][j][k] - f_next[i][j][k]);
                    f[i][j][k] = f_next[i][j][k];
                }
            }
        }
        residual = sum / (N*N*Q);

        iter++;
    }

    std::cout << "Simulation completed in " << iter << " iterations with a residual of " << residual << std::endl;

    return 0;
}