#include <iostream>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
#include<chrono>

__global__ void initialize(float* f, int* isn, float* wt, int nx, int ny, float rho0) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int a;

    if (idx < nx && idy < ny) {
        isn[idx * ny + idy] = 0;
        if (idy == 0 || idy == ny - 1) {
            isn[idx * ny + idy] = 1;
        }
        for (a = 0; a < 9; a++) {
            f[(a * nx * ny) + (idx * ny) + idy] = wt[a] * rho0;
        }
    }
}

__global__ void simulationStep(float* f, float* ft, int* isn, float* wt, int* ex, int* ey,
                              int* kb, float* source, float* feq, float dpdx, float tau,
                              int nx, int ny) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int a, ia, ja;
    float ux, uy, rho, term1, term2, u2;

    if (idx < nx && idy < ny) {
        ux = 0.0;
        uy = 0.0;
        rho = 0.0;
        for (a = 0; a < 9; a++) {
            rho += f[(a * nx * ny) + (idx * ny) + idy];
            ux += f[(a * nx * ny) + (idx * ny) + idy];
            uy += f[(a * nx * ny) + (idx * ny) + idy];
        }
        ux += dpdx / 2.0;
        ux /= rho;
        uy /= rho;
        u2 = ux * ux + uy * uy;

        for (a = 0; a < 9; a++) {
            term1 = ux * ex[a] + uy * ey[a];
            term2 = term1 * term1;
            source[a] = (1.0 - 0.5 / tau) * wt[a] * (3.0 * (ex[a] - ux) + 9.0 * (ex[a] * ux + ey[a] * uy) * ex[a]) * dpdx;
            feq[a] = wt[a] * rho * (1.0 + 3.0 * term1 + 4.5 * term2 - 1.5 * u2);
            ft[(a * nx * ny) + (idx * ny) + idy] = f[(a * nx * ny) + (idx * ny) + idy] - (f[(a * nx * ny) + (idx * ny) + idy] - feq[a]) / tau + source[a];
        }

        __syncthreads();

        ia = idx + ex[a];
        ja = idy + ey[a];
        if (ia < 0) {
            ia = nx - 1;
        }
        if (ia > nx - 1) {
            ia = 0;
        }
        f[(a * nx * ny) + (ia * ny) + ja] = ft[(a * nx * ny) + (idx * ny) + idy];
    }
}

int main() {
    //auto Tstart_before_KC = std::chrono::high_resolution_clock::now();
    int nx = 5001, ny = 2001;
    int a;
    float* f, * ft, * source, * feq;
    int* isn;
    float* wt;
    int* ex, * ey, * kb;
    float rho0 = 1.0;
    float dpdx = 1.0e-5;
    float tau = 0.8;
    int time = 2000;

    auto Tstart_malloc = std::chrono::high_resolution_clock::now();
// Allocate CPU memory
    f = new float[9 * nx * ny];
    ft = new float[9 * nx * ny];
    isn = new int[nx * ny];
    wt = new float[9];
    ex = new int[9];
    ey = new int[9];
    kb = new int[9];
    source = new float[9];
    feq = new float[9];

    // Allocate GPU memory
    float* d_f, * d_ft, * d_source, * d_feq;
    int* d_isn;
    float* d_wt;
    int* d_ex, * d_ey, * d_kb;

    cudaMalloc((void**)&d_f, 9 * nx * ny * sizeof(float));
    cudaMalloc((void**)&d_ft, 9 * nx * ny * sizeof(float));
    cudaMalloc((void**)&d_isn, nx * ny * sizeof(int));
    cudaMalloc((void**)&d_wt, 9 * sizeof(float));
    cudaMalloc((void**)&d_ex, 9 * sizeof(int));
    cudaMalloc((void**)&d_ey, 9 * sizeof(int));
    cudaMalloc((void**)&d_kb, 9 * sizeof(int));
    cudaMalloc((void**)&d_source, 9 * sizeof(float));
    cudaMalloc((void**)&d_feq, 9 * sizeof(float));
    auto Tend_malloc = std::chrono::high_resolution_clock::now();
    auto time_duration_malloc = std::chrono::duration_cast<std::chrono::milliseconds>(Tend_malloc - Tstart_malloc).count();
    std::cout << "Memory allocation time: " << time_duration_malloc << " milliseconds" << std::endl;

   // Initialize weights and velocities
    for (a = 0; a < 9; a++) {
        if (a == 0) {
            wt[a] = 4.0 / 9.0;
        }
        if (a >= 1 && a <= 4) {
            wt[a] = 1.0 / 9.0;
        }
        if (a >= 5 && a <= 8) {
            wt[a] = 1.0 / 36.0;
        }
    }

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

    for (int a = 0; a < 9; a++) {
        for (int a1 = a; a1 < 9; a1++) {
            if ((ex[a] + ex[a1] == 0 && ey[a] + ey[a1] == 0)) {
                kb[a] = a1;
                kb[a1] = a;
            }
        }
    }

    // Set grid and block dimensions for GPU kernels
    dim3 blockSize(16, 16);
    dim3 gridSize((nx + blockSize.x - 1) / blockSize.x, (ny + blockSize.y - 1) / blockSize.y);
//auto Tend_before_KC = std::chrono::high_resolution_clock::now();
//auto time_duration_before_KC = std::chrono::duration_cast<std::chrono::milliseconds>(Tend_before_KC - Tstart_before_KC).count();
//std::cout<< time_duration_before_KC << " milliseconds " << std::endl;_malloc

auto TstartSL = std::chrono::high_resolution_clock::now();
    for (int ts = 1; ts <= time; ts++) {
        auto TstartKernel = std::chrono::high_resolution_clock::now();
        simulationStep<<<gridSize, blockSize>>>(d_f, d_ft, d_isn, d_wt, d_ex, d_ey, d_kb, d_source, d_feq, dpdx, tau, nx, ny);
        cudaDeviceSynchronize();
        auto TendKernel = std::chrono::high_resolution_clock::now();
        auto time_duration_kernel = std::chrono::duration_cast<std::chrono::milliseconds>(TendKernel - TstartKernel).count();
        std::cout << "Simulation step " << ts << " kernel time: " << time_duration_kernel << " milliseconds" << std::endl;
    }
auto TendSL = std::chrono::high_resolution_clock::now();
auto time_durationSL = std::chrono::duration_cast<std::chrono::milliseconds>(TendSL - TstartSL).count();
std::cout << "Simulation loop time: " << time_durationSL << " milliseconds" << std::endl;

    // Copy data from GPU to CPU
    cudaMemcpy(f, d_f, 9 * nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(ft, d_ft, 9 * nx * ny * sizeof(float), cudaMemcpyDeviceToHost);
    delete[] f;
    delete[] ft;
    delete[] isn;
    delete[] wt;
    delete[] ex;
    delete[] ey;
    delete[] kb;
    delete[] source;
    delete[] feq;

    // Free GPU memory
    cudaFree(d_f);
    cudaFree(d_ft);
    cudaFree(d_isn);
    cudaFree(d_wt);
    cudaFree(d_ex);
    cudaFree(d_ey);
    cudaFree(d_kb);
    cudaFree(d_source);
    cudaFree(d_feq);

    return 0;
}