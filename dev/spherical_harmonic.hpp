#ifndef SPHERICAL_HARMONIC_HPP
#define SPHERICAL_HARMONIC_HPP

#include "include/odin/models/gravitational/gravitational_base.hpp"
#include "include/odin/models/gravitational/spherical.hpp"

#include <array>
#include <cassert>
#include <utility>
#include <vector>

template<typename T, std::size_t lmax = 0, std::size_t mmax = 0>
struct SphericalHarmonics {
    std::array<std::array<T, mmax + 1>, lmax + 1> Cnm{};
    std::array<std::array<T, mmax + 1>, lmax + 1> Snm{};

    // This constructor is only valid for statically sized harmonics
    SphericalHarmonics() {
        static_assert(lmax > 0 && mmax > 0, "lmax and mmax must be greater than 0 for statically sized harmonics");
    }

    // This constructor is for dynamically sized harmonics
    explicit SphericalHarmonics(
            std::size_t lsize, std::size_t msize) {
        assert(lsize > 0 && msize > 0);
        Cnm.assign(lsize + 1, std::vector<T>(msize + 1, 0));
        Snm.assign(lsize + 1, std::vector<T>(msize + 1, 0));
    }

    void set_coefficients(
            std::array<std::array<T, mmax + 1>, lmax + 1> Cnm_values,
            std::array<std::array<T, mmax + 1>, lmax + 1> Snm_values) {
        Cnm = std::move(Cnm_values);
        Snm = std::move(Snm_values);
    }

    void set_coefficients(std::vector<std::vector<T>> Cnm_values,
                          std::vector<std::vector<T>> Snm_values) {
        Cnm = std::move(Cnm_values);
        Snm = std::move(Snm_values);
    }
};

template<typename S, size_t Degree, size_t Order, size_t Dim = 3>
class [[maybe_unused]] SphericalHarmonicBase : public SphericalBase<SphericalHarmonicBase<S, Degree, Order, Dim>, S, Dim> {
    using Base = SphericalBase<SphericalHarmonicBase<S, Degree, Order, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit SphericalHarmonicBase(
            const Scalar                                     mu,
            const int                                        max_degree,
            const Scalar                                     reference_radius,
            const SphericalHarmonics<Scalar, Degree, Order> &spherical_harmonics,
            const Vector                                    &position = Vector::Zero())
        : Base(mu, reference_radius, position) {

        spherical_harmonics_ = std::make_unique<SphericalHarmonics<Scalar, Degree, Order>>(spherical_harmonics);
    }

    Scalar static_potential_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        Scalar r     = rel_position.norm();
        Scalar theta = acos(rel_position[2] / r);              // z/r
        Scalar phi   = atan2(rel_position[1], rel_position[0]);// y/x
        return compute_gravitational_potential(r, theta, phi, mu, radius);
    }

    Vector static_acceleration_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        Scalar r     = rel_position.norm();
        Scalar theta = acos(rel_position[2] / r);              // z/r
        Scalar phi   = atan2(rel_position[1], rel_position[0]);// y/x
        return spherical_to_cartesian(compute_acceleration(r, theta, phi, mu, radius), theta, phi);
    }

private:
    shtns_cfg                                                  shtns_;
    std::unique_ptr<SphericalHarmonics<Scalar, Degree, Order>> spherical_harmonics_{};

    Vector spherical_to_cartesian(const Vector &spherical, Scalar theta, Scalar phi) {
        // Convert from spherical coordinates (r, theta, phi) to cartesian coordinates (x, y, z)
        Scalar sin_theta = sin(theta);
        Scalar cos_theta = cos(theta);
        Vector cartesian;
        cartesian[0] = spherical[0] * sin_theta * cos(phi) + spherical[1] * cos_theta * cos(phi) - spherical[2] * sin(phi);// x-component
        cartesian[1] = spherical[0] * sin_theta * sin(phi) + spherical[1] * cos_theta * sin(phi) + spherical[2] * cos(phi);// y-component
        cartesian[2] = spherical[0] * cos_theta - spherical[1] * sin_theta;                                                // z-component
        return cartesian;
    }


    Scalar compute_gravitational_potential(Scalar r, Scalar theta, Scalar phi, Scalar mu, Scalar radius) {
        Scalar U = 0;
        //        int    lmax = spherical_harmonics_->get_lmax();
        //        int    mmax = spherical_harmonics_->get_mmax();
        for (int n = 0; n <= lmax_; ++n) {
            for (int m = 0; m <= std::min(n, mmax_); ++m) {
                Scalar Pnm    = shtns_legendre_nm(n, m, std::cos(theta));// Using cos(theta) because theta is co-latitude
                Scalar factor = std::pow(radius / r, n + 1) * Pnm;
                U += factor * (spherical_harmonics_->Cnm(n, m) * std::cos(m * phi) + spherical_harmonics_->Snm(n, m) * std::sin(m * phi));
            }
        }
        return mu / r * U;
    }

    Vector compute_acceleration(Scalar r, Scalar theta, Scalar phi, Scalar mu, Scalar radius) {
        Vector acc = Vector::Zero();
        // lmax - maximum degree of spherical harmonics
        // mmax - maximum order of spherical harmonics
        for (int n = 0; n <= lmax; ++n) {
            for (int m = 0; m <= std::min(n, mmax); ++m) {
                Scalar Pnm      = std::assoc_legendre(n, m, std::cos(theta));// Legendre function
                Scalar dPnm     = std::assoc_legendre(n, m, std::cos(theta));// Derivative of Legendre function
                Scalar factor   = std::pow(radius / r, n + 1);
                Scalar cos_mphi = std::cos(m * phi);
                Scalar sin_mphi = std::sin(m * phi);
                acc[0] += factor * (n + 1) * Pnm * (spherical_harmonics_->Cnm(n, m) * cos_mphi + spherical_harmonics_->Snm(n, m) * sin_mphi);
                acc[1] += factor * (dPnm / std::sin(theta)) * (spherical_harmonics_->Cnm(n, m) * cos_mphi + spherical_harmonics_->Snm(n, m) * sin_mphi);
                acc[2] += factor * m * Pnm * (-spherical_harmonics_->Cnm(n, m) * sin_mphi + spherical_harmonics_->Snm(n, m) * cos_mphi);
            }
        }
        acc[0] = -mu / (r * r) * acc[0];
        acc[1] = -mu / r * acc[1];
        acc[2] = -mu / (r * std::sin(theta)) * acc[2];
        return acc;
    }
};

#endif// SPHERICAL_HARMONIC_HPP