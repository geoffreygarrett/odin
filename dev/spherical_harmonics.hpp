#include <cmath>
#include <shtns.h>
#include <stdio.h>
#include <vector>

#define GM 3.986004418e14// Earth's gravitational constant in m^3/s^2
#define R 6.3781e6       // Earth's radius in m
#define MAX_DEGREE 8     // maximum degree for spherical harmonics

// precomputed constants for the potential model
std::vector<std::vector<double>> Cnm, Snm;

// your SHTns configuration
shtns_cfg shtns;

//int num_threads = 4;  // or however many you want to use
//shtns_use_threads(num_threads);
void initialize_spherical_harmonic_coefficients(int lmax, int mmax) {
    // Initialize all coefficients to zero
    Cnm.assign(lmax + 1, std::vector<double>(mmax + 1, 0));
    Snm.assign(lmax + 1, std::vector<double>(mmax + 1, 0));

    // Zonal coefficients
    Cnm[2][0] = -0.1082635854e-2;
    Cnm[3][0] = 0.2532435346e-5;
    Cnm[4][0] = 0.1619331205e-5;
    Cnm[5][0] = 0.2277161016e-6;
    Cnm[6][0] = -0.5396484906e-6;
    Cnm[7][0] = 0.3513684422e-6;
    Cnm[8][0] = 0.2025187152e-6;

    // Tesseral coefficients
    Cnm[2][1] = -0.3504890360e-9;
    Cnm[2][2] = 0.1574536043e-5;
    Cnm[3][1] = 0.2192798802e-5;
    Cnm[3][2] = 0.3090160446e-6;
    Cnm[3][3] = 0.1005588574e-6;
    Cnm[4][1] = -0.5087253036e-6;
    Cnm[4][2] = 0.7841223074e-7;
    Cnm[4][3] = 0.5921574319e-7;
    Cnm[4][4] = -0.3982395740e-8;

    Snm[2][1] = 0.1635406077e-8;
    Snm[2][2] = -0.9038680729e-6;
    Snm[3][1] = 0.2680118938e-6;
    Snm[3][2] = -0.2114023978e-6;
    Snm[3][3] = 0.1972013239e-6;
    Snm[4][1] = -0.4494599352e-6;
    Snm[4][2] = 0.1481554569e-6;
    Snm[4][3] = -0.1201129183e-7;
    Snm[4][4] = 0.6525605810e-8;
}


double compute_gravitational_potential(double r, double theta, double phi) {
    double U    = 0;
    int    lmax = Cnm.size() - 1;
    int    mmax = Cnm[0].size() - 1;
    for (int n = 0; n <= lmax; ++n) {
        for (int m = 0; m <= std::min(n, mmax); ++m) {
            double Pnm    = shtns_legendre_nm(n, m, std::cos(theta));// Using cos(theta) because theta is co-latitude
            double factor = std::pow(R / r, n + 1) * Pnm;
            U += factor * (Cnm[n][m] * std::cos(m * phi) + Snm[n][m] * std::sin(m * phi));
        }
    }
    return GM / r * U;
}

template<typename T>
T compute_gravitational_potential(T r, T theta, T phi) {
    T   U    = 0;
    int lmax = Cnm.size() - 1;
    int mmax = Cnm[0].size() - 1;

    tbb::parallel_for(0, lmax + 1, [&](int n) {
        for (int m = 0; m <= std::min(n, mmax); ++m) {
            T Pnm    = shtns_legendre_nm<T>(n, m, std::cos(theta));// Assuming these functions are template functions.
            T factor = std::pow(R / r, n + 1) * Pnm;
            tbb::atomic_fetch_and_add(&U, factor * (Cnm[n][m] * std::cos(m * phi) + Snm[n][m] * std::sin(m * phi)));
        }
    });

    return GM / r * U;
}


typename<typename T, bool use_tbb = false>
        std::vector<double>

        compute_acceleration(double r, double theta, double phi) {
    double Ur = 0, Utheta = 0, Uphi = 0;
    int    lmax = Cnm.size() - 1;
    int    mmax = Cnm[0].size() - 1;
    for (int n = 0; n <= lmax; ++n) {
        for (int m = 0; m <= std::min(n, mmax); ++m) {
            double Pnm      = shtns_legendre_nm(n, m, std::cos(theta)); // Legendre function
            double dPnm     = shtns_dlegendre_nm(n, m, std::cos(theta));// Derivative of Legendre function
            double factor   = std::pow(R / r, n + 1);
            double cos_mphi = std::cos(m * phi);
            double sin_mphi = std::sin(m * phi);
            Ur += factor * (n + 1) * Pnm * (Cnm[n][m] * cos_mphi + Snm[n][m] * sin_mphi);
            Utheta += factor * (dPnm / std::sin(theta)) * (Cnm[n][m] * cos_mphi + Snm[n][m] * sin_mphi);
            Uphi += factor * m * Pnm * (-Cnm[n][m] * sin_mphi + Snm[n][m] * cos_mphi);
        }
    }
    double a_r     = -GM / (r * r) * Ur;
    double a_theta = -GM / r * Utheta;
    double a_phi   = -GM / (r * std::sin(theta)) * Uphi;
    return {a_r, a_theta, a_phi};
}

template<typename T, bool use_tbb = true>
std::vector<T> compute_acceleration(T r, T theta, T phi) {
    T   Ur = 0, Utheta = 0, Uphi = 0;
    int lmax = Cnm.size() - 1;
    int mmax = Cnm[0].size() - 1;

    tbb::parallel_for(0, lmax + 1, [&](int n) {
        for (int m = 0; m <= std::min(n, mmax); ++m) {
            T Pnm      = shtns_legendre_nm<T>(n, m, std::cos(theta));// Assuming these functions are template functions.
            T dPnm     = shtns_dlegendre_nm<T>(n, m, std::cos(theta));
            T factor   = std::pow(R / r, n + 1);
            T cos_mphi = std::cos(m * phi);
            T sin_mphi = std::sin(m * phi);
            tbb::atomic_fetch_and_add(
                    &Ur,
                    factor * (n + 1) * Pnm * (Cnm[n][m] * cos_mphi + Snm[n][m] * sin_mphi));
            tbb::atomic_fetch_and_add(
                    &Utheta,
                    factor * (dPnm / std::sin(theta)) * (Cnm[n][m] * cos_mphi + Snm[n][m] * sin_mphi));
            tbb::atomic_fetch_and_add(
                    &Uphi,
                    factor * m * Pnm * (-Cnm[n][m] * sin_mphi + Snm[n][m] * cos_mphi));
        }
    });

    T a_r     = -GM / (r * r) * Ur;
    T a_theta = -GM / r * Utheta;
    T a_phi   = -GM / (r * std::sin(theta)) * Uphi;
    return {a_r, a_theta, a_phi};
}

int main() {
    int lmax = MAX_DEGREE;// maximum degree
    int mmax = lmax;      // maximum order

    // Initialize the SHTns configuration
    shtns = shtns_init(sht_reg_fast, lmax, mmax, 0, sht_orthonormal + sht_quick_init);

    // Initialize spherical harmonic coefficients
    initialize_spherical_harmonic_coefficients(lmax, mmax);

    // Compute the gravitational potential and acceleration at some point
    double              r            = R + 1000;// altitude of 1km
    double              theta        = M_PI / 4;// 45 degrees latitude
    double              phi          = M_PI / 4;// 45 degrees longitude
    double              U            = compute_gravitational_potential(r, theta, phi);
    std::vector<double> acceleration = compute_acceleration(r, theta, phi);

    // Do something with U and acceleration
    printf("Potential: %f\n", U);
    printf("Acceleration: %f, %f, %f\n", acceleration[0], acceleration[1], acceleration[2]);

    return 0;
}