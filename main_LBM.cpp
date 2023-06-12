#include<iostream>
#include<cmath>

const int grid_size = 50;
const double L = 1.0;
const double max_vel = 0.1;
const double nu = 0.1;
const double Rho_0 = 1.0;
const double dt = 1.0;
const double dx = L/grid_size;
const int maximum_iterations = 1000;
const double converge = 1e-5;

//Lattice properties
const int D2Q9 = 9;
const int C_x[D2Q9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int C_y[D2Q9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double Weights[D2Q9] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };
double Tau = 3.0*nu + 0.5;
double omega = 1.0/Tau;

//Distribution function Declaration and Initialization
double f[grid_size][grid_size][D2Q9];
double f_new[grid_size][grid_size][D2Q9];

void initialize_f(){
    for(int i=0; i< grid_size; ++i){
        for(int j=0; j< grid_size; ++j){
            double u = max_vel * 6.0 * (j*(grid_size - j))/ (grid_size * grid_size);
            double Rho = Rho_0;

            for(int k=0; k< D2Q9; ++k){
                double C_u = 3.0 * (C_x[k] * u);
                double U_u = 3.0 * (C_x[k] * u) * (C_x[k] * u);
                f[i][j][k] = Weights[k] * Rho * (1.0 + C_u + 0.5*U_u - 1.5*u*u);
            }
        }
    }
}


//Collision and Streaming Process Function
void collision_stream(){
    //1. Collision
    for(int i=0; i< grid_size; ++i){
        for(int j=0; j< grid_size; ++j){
            double Rho = 0.0, U_x = 0.0, U_y = 0.0;

            for(int k= 0; k< D2Q9; ++k){
                Rho += f[i][j][k];
                U_x += C_x[k]*f[i][j][k];
                U_y += C_y[k]*f[i][j][k];
            }

            U_x /= Rho;
            U_y /= Rho;

            for(int k=0; k< D2Q9; ++k){
                double C_u = 3.0 * (C_x[k]*U_x + C_y[k]*U_y);
                double U_u = 3.0 * (C_x[k] * U_x + C_y[k] * U_y) * (C_x[k] * U_x + C_y[k] * U_y);
                f_new[i][j][k] = f [i][j][k] - (f[i][j][k] - Weights[k]*Rho*(1.0 + C_u + 0.5 * U_u - 1.5 * (U_x * U_x + U_y * U_y))) / Tau; 
            }
        }
    }

    //2. Streaming
        for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            for (int k = 0; k < D2Q9; ++k) {
                int i_new = (i + C_x[k] + grid_size) % grid_size;
                int j_new = (j + C_y[k] + grid_size) % grid_size;
                f[i_new][j_new][k] = f_new[i][j][k];
            }
        }
    }
}

//Calculation of Maximum Velocity
double Max_Velocity(){
    double maxVel = 0.0;
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            double U_x = 0.0, U_y = 0.0;

            for (int k = 0; k < D2Q9; ++k) {
                U_x += C_x[k] * f[i][j][k];
                U_y += C_y[k] * f[i][j][k];
            }

            double vel = std::sqrt(U_x * U_x + U_y * U_y);
            if (vel > maxVel) {
                maxVel = vel;
            }
        } 
    }
    return maxVel;
}


int main(){
    initialize_f();
    int counter = 0;
    double residual = converge + 1.0;

    for(int counter= 0; counter< maximum_iterations; counter++){
        double MaxVelOut = Max_Velocity();
        double U_conv = max_vel* dx / MaxVelOut;
        const int dt = dx*dx/(U_conv * U_conv);

        collision_stream();
        //std::cout<< MaxVelOut << std::endl;
        double sum=0.0;
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                for (int k = 0; k < D2Q9; ++k) {
                    sum += std::abs(f[i][j][k] - f_new[i][j][k]);
                    f[i][j][k] = f_new[i][j][k];
                }
            }
        }
    residual = sum/(grid_size*grid_size*D2Q9);

    }

/*
    while(counter < maximum_iterations && residual > converge){
        double MaxVelOut = Max_Velocity();
        double U_conv = max_vel* dx / MaxVelOut;
        const int dt = dx*dx/(U_conv * U_conv);

        collision_stream();
        //std::cout<< MaxVelOut << std::endl;
        double sum=0.0;
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                for (int k = 0; k < D2Q9; ++k) {
                    sum += std::abs(f[i][j][k] - f_new[i][j][k]);
                    f[i][j][k] = f_new[i][j][k];
                }
            }
        }
    residual = sum/(grid_size*grid_size*D2Q9);
    counter++;
    }
*/
   
    return 0;
}

