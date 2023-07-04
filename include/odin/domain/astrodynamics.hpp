#ifndef ASTRODYNAMICS_HPP
#define ASTRODYNAMICS_HPP

#include <Eigen/Dense>


#include <cmath>
//#include <ranges>
#include <type_traits>

#define ODIN_USE_TBB
#ifdef ODIN_USE_TBB

#include <tbb/parallel_for.h>

#endif


template<typename Float>
Float PI = static_cast<Float>(3.14159265358979323846264338327950288419716939937510L);


template<typename Float>
using Vector3 = Eigen::Vector3<Float>;

template<typename T>
T mod(T a, T b) {
    return std::fmod((std::fmod(a, b) + b), b);
}

template<typename Float>
Float anomaly_eccentric_to_mean(Float E, Float e) {
    return E - e * std::sin(E);
}

template<typename Float>
Float anomaly_true_to_eccentric(Float f, Float e) {
    return Float(2) * std::atan(std::sqrt((Float(1) - e) / (Float(1) + e)) * std::tan(f / Float(2)));
}

template<typename Float>
Float anomaly_eccentric_to_true(Float E, Float e) {
    Float beta = e / (Float(1) + std::sqrt(Float(1) - e * e));
    return E + Float(2) * std::atan(beta * std::sin(E) / (Float(1) - beta * std::cos(E)));
}

template<typename Float>
Float anomaly_mean_to_eccentric(Float M, Float e, int max_iter = 1000, Float tolerance = 1e-6) {
    // improved initial guess
    Float E = (e < Float(0.8)) ? M : PI<Float>;
    Float delta = std::numeric_limits<Float>::max();
    int iter = 0;

    // Newton-Raphson method to solve Kepler's equation
    while (delta > tolerance && iter < max_iter) {
        Float E_next = E - (E - e * std::sin(E) - M) / (Float(1.) - e * std::cos(E));
        delta = std::abs(E_next - E);
        E = E_next;
        ++iter;
    }

    return E;
}

template<typename Float>
Float anomaly_mean_to_true(Float M, Float e, int max_iter = 1000, Float tolerance = 1e-6) {
    // Convert from mean anomaly to eccentric anomaly
    Float E = anomaly_mean_to_eccentric(M, e, max_iter, tolerance);

    // Convert from eccentric anomaly to true anomaly
    return anomaly_eccentric_to_true(E, e);
}

template<typename Float>
std::vector<Float>
anomaly_mean_to_true(std::vector<Float> MM, std::vector<Float> ee, int max_iter = 1000, Float tolerance = 1e-6) {
    // Convert from mean anomaly to eccentric anomaly
    std::vector<Float> EE(MM.size());
    for (int i = 0; i < MM.size(); ++i) {
        EE[i] = anomaly_mean_to_eccentric(MM[i], ee[i], max_iter, tolerance);
    }

    // Convert from eccentric anomaly to true anomaly
    std::vector<Float> ff(MM.size());
    for (int i = 0; i < MM.size(); ++i) {
        ff[i] = anomaly_eccentric_to_true(EE[i], ee[i]);
    }

    return ff;
}


template<typename Float>
Float anomaly_true_to_mean(Float f, Float e) {
    // Convert from true anomaly to eccentric anomaly
    Float E = anomaly_true_to_eccentric(f, e);

    // Convert from eccentric anomaly to mean anomaly
    return anomaly_eccentric_to_mean(E, e);
}

template<typename Float>
Float anomaly_hyperbolic_to_mean(Float F, Float e) {
    return e * std::sinh(F) - F;
}

template<typename Float>
// TODO: Check this
Float anomaly_mean_to_hyperbolic(Float M, Float e, int max_iter = 1000, Float tolerance = 1e-6) {
    Float F = M;
    Float dF = Float(0);
    for (int i = 0; i < max_iter; ++i) {
        dF = (e * std::sinh(F) - F - M) / (e * std::cosh(F) - Float(1));
        F -= dF;
        if (std::abs(dF) < tolerance) {
            break;
        }
    }
    return F;
}

template<typename Float>
// TODO: Check this
Float anomaly_true_to_hyperbolic(Float f, Float e) {
    return Float(2) * std::atan(std::sqrt((e - Float(1)) / (e + Float(1))) * std::tan(f / Float(2)));
}

template<typename Float>
// TODO: Check this
Float anomaly_hyperbolic_to_true(Float F, Float e) {
    Float beta = e / (Float(1) + std::sqrt(e * e - Float(1)));
    return Float(2) * std::atan(beta * std::sinh(F / Float(2)) / (Float(1) + beta * std::cosh(F / Float(2))));
}

/**
 * @brief Converts Cartesian to Classical Orbital Elements (COE)
 *
 * @param r Position vector in Cartesian coordinates
 * @param v Velocity vector in Cartesian coordinates
 * @param mu Gravitational constant
 * @return std::tuple<Float, Float, Float, Float, Float, Float> Classical Orbital Elements (p, a, e, i, Omega, nu)
 *
 * @details This algorithm is based on the RV2COE algorithm presented in:
 * David A. Vallado, James Wertz - Fundamentals of Astrodynamics and Applications, 4th ed (Page 114)
 *
 */
template<typename Float>
std::tuple<Float, Float, Float, Float, Float, Float> rv2coe(
        Float mu,
        const Vector3<Float> &r,
        const Vector3<Float> &v,
        Float tol = 1e-8) {
    // h = r x v, h_mag = |h|
    Vector3<Float> h = r.cross(v);
    Float h_norm = h.norm();

    // n = k_hat x h, n_mag = |n|
    Vector3<Float> n = Vector3<Float>::UnitZ().cross(h);
    Float n_norm = n.norm();

    // e = ((v^2 - mu/r)*r - (r.v)*v)/mu, e_mag = |e|
    Vector3<Float> e = (r * (v.squaredNorm() - mu / r.norm()) - v * r.dot(v)) / mu;
    Float e_norm = e.norm();

    // p = h^2 / mu
    Float p = h.squaredNorm() / mu;

    // i = acos(h_z / |h|)
    Float inc = std::acos(h.z() / h_norm);

    // Special Cases
    bool circular = e_norm < tol;
    bool equatorial = std::abs(inc) < tol;


    Float raan, argp, nu;

    if (equatorial && !circular) {
        raan = 0;
        argp = mod(std::atan2(e.y(), e.x()), 2 * PI<Float>);
        nu = std::atan2(h.dot(e.cross(r)) / h_norm, e.dot(r));
    } else if (!equatorial && circular) {
        raan = mod(std::atan2(n.y(), n.x()), 2 * PI<Float>);
        argp = 0;
        nu = std::atan2(r.dot(h.cross(n)) / h_norm, r.dot(n));
    } else if (equatorial && circular) {
        raan = 0;
        argp = 0;
        nu = mod(std::atan2(r.y(), r.x()), 2 * PI<Float>);
    } else {
        Float a = p / (1 - std::pow(e_norm, 2));
        Float mu_a = mu * a;
        if (a > 0) {
            Float e_se = r.dot(v) / std::sqrt(mu_a);
            Float e_ce = r.norm() * v.squaredNorm() / mu - 1;
            nu = anomaly_eccentric_to_true(std::atan2(e_se, e_ce), e_norm);
        } else {
            Float e_sh = r.dot(v) / std::sqrt(-mu_a);
            Float e_ch = r.norm() * v.squaredNorm() / mu - 1;
            nu = anomaly_hyperbolic_to_true(std::log((e_ch + e_sh) / (e_ch - e_sh)) / 2, e_norm);
        }
        raan = mod(std::atan2(n.y(), n.x()), 2 * PI<Float>);
        Float px = r.dot(n);
        Float py = r.dot(h.cross(n)) / h_norm;
        argp = mod(std::atan2(py, px) - nu, 2 * PI<Float>);
    }
    nu = mod(nu + PI<Float>, 2 * PI<Float>) - PI<Float>;
    return std::make_tuple(p, e_norm, inc, raan, argp, nu);
}

template<typename Float>
constexpr auto rv2coe(
        Float mu,
        const std::vector<Vector3<Float>> &rr,
        const std::vector<Vector3<Float>> &vv,
        Float tol = 1e-8) {
    std::vector<Float> p(rr.size());
    std::vector<Float> e(rr.size());
    std::vector<Float> inc(rr.size());
    std::vector<Float> raan(rr.size());
    std::vector<Float> argp(rr.size());
    std::vector<Float> nu(rr.size());
    for (size_t i = 0; i < rr.size(); ++i) {
        std::tie(p[i], e[i], inc[i], raan[i], argp[i], nu[i]) = rv2coe(mu, rr[i], vv[i], tol);
    }
    return std::make_tuple(p, e, inc, raan, argp, nu);
}

template<typename Float>
constexpr std::tuple<Eigen::Vector3 < Float>, Eigen::Vector3 <Float>>
coe2pqr(
        Float
p,
Float e,
        Float
nu,
Float mu
) {
Float cos_nu = std::cos(nu);
Float sin_nu = std::sin(nu);
Float sqrt_mu_p = std::sqrt(mu / p);
return
std::make_tuple(
        Eigen::Vector3<Float>{
                p * cos_nu / (1 + e * cos_nu),
                p * sin_nu / (1 + e * cos_nu),
                0},
        Eigen::Vector3<Float>{
        -sqrt_mu_p * sin_nu,
        sqrt_mu_p * (e + cos_nu),
        0}
);
}

template<int N, typename Float>
constexpr Eigen::Matrix3 <Float> rot(const Float &angle) {
    static_assert(N >= 0 && N <= 2, "N must be between 0 and 2 (inclusive and 0-indexed)");
    const Eigen::Vector3 <Float> unit_vectors[] = {
            Eigen::Vector3<Float>::UnitX(),
            Eigen::Vector3<Float>::UnitY(),
            Eigen::Vector3<Float>::UnitZ()};
    return Eigen::AngleAxis<Float>(angle, unit_vectors[N]).toRotationMatrix();
}

template<typename Float>
constexpr std::tuple<Eigen::Vector3 < Float>, Eigen::Vector3 <Float>>
coe2rv(
        Float
mu,  // gravitational parameter
Float p,   // semi-latus rectum
Float
e,   // eccentricity
Float i,   // inclination
Float
raan,// right ascension of ascending node
Float argp,// argument of periapsis
Float
nu   // true anomaly
) {
// Need to implement the coe2pqr function, not provided in the original code
auto [r_pqw, v_pqw] = coe2pqr<Float>(p, e, nu, mu);

// create Eigen rotation matrices
Eigen::Matrix3 <Float> rot_ijk_pqw = rot<2, Float>(raan)
                                     * rot<0, Float>(i)
                                     * rot<2, Float>(argp);

// return the rotated vectors
return
std::make_tuple(
        rot_ijk_pqw
* r_pqw,
rot_ijk_pqw *v_pqw
);
}

template<typename Float>
constexpr std::tuple<Eigen::Matrix < Float, 3, Eigen::Dynamic>, Eigen::Matrix<Float, 3, Eigen::Dynamic>>

coe2rv(
        const Float mu,  // gravitational parameter
        const Float p,   // semi-latus rectum
        const Float e,   // eccentricity
        const Float inc, // inclination
        const Float raan,// right ascension of ascending node
        const Float argp,// argument of periapsis
        const std::vector<Float> &nu   // true anomaly
) {
    size_t n = nu.size();
    Eigen::Matrix<Float, 3, Eigen::Dynamic> r_out(3, n);
    Eigen::Matrix<Float, 3, Eigen::Dynamic> v_out(3, n);

#ifdef ODIN_USE_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range <size_t> &range) {
        for (size_t i = range.begin(); i != range.end(); ++i) {
            auto temp = coe2rv(mu, p, e, inc, raan, argp, nu[i]);
            r_out.col(i) = std::get<0>(temp);
            v_out.col(i) = std::get<1>(temp);
        }
    });
#else
    for (size_t i = 0; i < n; ++i) {
        auto temp    = coe2rv(mu, p, e, inc, raan, argp, nu[i]);
        r_out.col(i) = std::get<0>(temp);
        v_out.col(i) = std::get<1>(temp);
    }
#endif
    return std::make_tuple(r_out, v_out);
}

template<typename Float>
constexpr std::tuple<Eigen::Matrix < Float, 3, Eigen::Dynamic>, Eigen::Matrix<Float, 3, Eigen::Dynamic>>

coe2rv(
        const Float mu,  // gravitational parameter
        const std::vector<Float> &p,   // semi-latus rectum
        const std::vector<Float> &e,   // eccentricity
        const std::vector<Float> &i,   // inclination
        const std::vector<Float> &raan,// right ascension of ascending node
        const std::vector<Float> &argp,// argument of periapsis
        const std::vector<Float> &nu   // true anomaly
) {
    size_t n = p.size();
    Eigen::Matrix<Float, 3, Eigen::Dynamic> r_out(3, n);
    Eigen::Matrix<Float, 3, Eigen::Dynamic> v_out(3, n);

#ifdef ODIN_USE_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range <size_t> &r) {
        for (size_t idx = r.begin(); idx != r.end(); ++idx) {
            auto [r_pqw, v_pqw] = coe2pqr<Float>(p[idx], e[idx], nu[idx], mu);
            Eigen::Matrix3 <Float> rot_ijk_pqw = rot<2, Float>(raan[idx])
                                                 * rot<0, Float>(i[idx])
                                                 * rot<2, Float>(argp[idx]);
            r_out.col(idx) = rot_ijk_pqw * r_pqw;
            v_out.col(idx) = rot_ijk_pqw * v_pqw;
        }
    });
#else
    for (size_t idx = 0; idx < n; ++idx) {
        auto [r_pqw, v_pqw]               = coe2pqr<Float>(p[idx], e[idx], nu[idx], mu);
        Eigen::Matrix3<Float> rot_ijk_pqw = rot<2, Float>(raan[idx])
                                          * rot<0, Float>(i[idx])
                                          * rot<2, Float>(argp[idx]);
        r_out.col(idx) = rot_ijk_pqw * r_pqw;
        v_out.col(idx) = rot_ijk_pqw * v_pqw;
    }
#endif
    return std::make_tuple(r_out, v_out);
}

template<typename Float, typename ConversionFunc>
auto sample_anomalies(Float e, int n_points, ConversionFunc conversion) {
    Float step = 2 * PI<Float>; / (n_points - 1);

    // Create a range, perform conversion
    // TODO: See if this is actually possible with LLVM 15/16. It's a bit of a hassle right now.
//    auto anomalies = std::views::iota(0, n_points)
//                   | std::views::transform([=](int i) { return i * step; })// Anomalies
//                   | std::views::transform(conversion);                    // Perform conversion
    std::vector<Float> anomalies(n_points);
    for (int i = 0; i < n_points; i++) {
        Float anomaly = i * step;  // Anomalies
        anomalies[i] = conversion(anomaly);  // Perform conversion
    }
    return std::vector<Float>(anomalies.begin(), anomalies.end());
}

template<typename Float>
auto sample_true_from_mean_anomaly(Float e, int n_points) {
    // Define conversion from mean to true anomaly
    auto conversion = [=](Float m) {
        return anomaly_eccentric_to_true(anomaly_mean_to_eccentric(m, e), e);
    };
    return sample_anomalies(e, n_points, conversion);
}

template<typename Float>
auto sample_true_from_eccentric_anomaly(Float e, int n_points) {
    // Define conversion from eccentric to true anomaly
    auto conversion = [=](Float E) {
        return anomaly_eccentric_to_true(E, e);
    };
    return sample_anomalies(e, n_points, conversion);
}

#endif// ASTRODYNAMICS_HPP
