#ifndef ASTRODYNAMICS_HPP
#define ASTRODYNAMICS_HPP

#include <Eigen/Dense>
#include <cmath>
//#include <ranges>
#include <type_traits>

#define ODIN_USE_TBB
#ifdef ODIN_USE_TBB

#include <odin/constants.h>
#include <tbb/parallel_for.h>
#endif

template<typename Scalar>
using Vector3 = Eigen::Vector3<Scalar>;

template<typename T>
T mod(T a, T b) {
    return std::fmod((std::fmod(a, b) + b), b);
}

namespace odin::domain::astrodynamics {

template<typename Scalar>
Scalar anomaly_eccentric_to_mean(Scalar E, Scalar e) {
    return E - e * std::sin(E);
}

template<typename Scalar>
Scalar anomaly_true_to_eccentric(Scalar f, Scalar e) {
    return Scalar(2)
         * std::atan(std::sqrt((Scalar(1) - e) / (Scalar(1) + e)) * std::tan(f / Scalar(2)));
}

template<typename Scalar>
Scalar anomaly_eccentric_to_true(Scalar E, Scalar e) {
    Scalar beta = e / (Scalar(1) + std::sqrt(Scalar(1) - e * e));
    return E + Scalar(2) * std::atan(beta * std::sin(E) / (Scalar(1) - beta * std::cos(E)));
}

template<typename Scalar>
Scalar anomaly_mean_to_eccentric(Scalar M, Scalar e, int max_iter = 1000, Scalar tolerance = 1e-6) {
    // improved initial guess
    Scalar E     = (e < Scalar(0.8)) ? M : odin::consts::PI<Scalar>;
    Scalar delta = std::numeric_limits<Scalar>::max();
    int    iter  = 0;

    // Newton-Raphson method to solve Kepler's equation
    while (delta > tolerance && iter < max_iter) {
        Scalar E_next = E - (E - e * std::sin(E) - M) / (Scalar(1.) - e * std::cos(E));
        delta         = std::abs(E_next - E);
        E             = E_next;
        ++iter;
    }

    return E;
}

template<typename Scalar>
Scalar anomaly_mean_to_true(Scalar M, Scalar e, int max_iter = 1000, Scalar tolerance = 1e-6) {
    // Convert from mean anomaly to eccentric anomaly
    Scalar E = anomaly_mean_to_eccentric(M, e, max_iter, tolerance);

    // Convert from eccentric anomaly to true anomaly
    return anomaly_eccentric_to_true(E, e);
}

template<typename Scalar>
std::vector<Scalar> anomaly_mean_to_true(std::vector<Scalar> MM,
        std::vector<Scalar>                                  ee,
        int                                                  max_iter  = 1000,
        Scalar                                               tolerance = 1e-6) {
    // Convert from mean anomaly to eccentric anomaly
    std::vector<Scalar> EE(MM.size());
    for (int i = 0; i < MM.size(); ++i) {
        EE[i] = anomaly_mean_to_eccentric(MM[i], ee[i], max_iter, tolerance);
    }

    // Convert from eccentric anomaly to true anomaly
    std::vector<Scalar> ff(MM.size());
    for (int i = 0; i < MM.size(); ++i) { ff[i] = anomaly_eccentric_to_true(EE[i], ee[i]); }

    return ff;
}


template<typename Scalar>
Scalar anomaly_true_to_mean(Scalar f, Scalar e) {
    // Convert from true anomaly to eccentric anomaly
    Scalar E = anomaly_true_to_eccentric(f, e);

    // Convert from eccentric anomaly to mean anomaly
    return anomaly_eccentric_to_mean(E, e);
}

template<typename Scalar>
Scalar anomaly_hyperbolic_to_mean(Scalar F, Scalar e) {
    return e * std::sinh(F) - F;
}

template<typename Scalar>
// TODO: Check this
Scalar anomaly_mean_to_hyperbolic(
        Scalar M, Scalar e, int max_iter = 1000, Scalar tolerance = 1e-6) {
    Scalar F  = M;
    Scalar dF = Scalar(0);
    for (int i = 0; i < max_iter; ++i) {
        dF = (e * std::sinh(F) - F - M) / (e * std::cosh(F) - Scalar(1));
        F -= dF;
        if (std::abs(dF) < tolerance) { break; }
    }
    return F;
}

template<typename Scalar>
// TODO: Check this
Scalar anomaly_true_to_hyperbolic(Scalar f, Scalar e) {
    return Scalar(2)
         * std::atan(std::sqrt((e - Scalar(1)) / (e + Scalar(1))) * std::tan(f / Scalar(2)));
}

template<typename Scalar>
// TODO: Check this
Scalar anomaly_hyperbolic_to_true(Scalar F, Scalar e) {
    Scalar beta = e / (Scalar(1) + std::sqrt(e * e - Scalar(1)));
    return Scalar(2)
         * std::atan(
                 beta * std::sinh(F / Scalar(2)) / (Scalar(1) + beta * std::cosh(F / Scalar(2))));
}

/**
 * @brief Converts Cartesian to Classical Orbital Elements (COE)
 *
 * @param r Position vector in Cartesian coordinates
 * @param v Velocity vector in Cartesian coordinates
 * @param mu Gravitational constant
 * @return std::tuple<Scalar, Scalar, Scalar, Scalar, Scalar, Scalar> Classical Orbital Elements (p, a, e, i, Omega, nu)
 *
 * @details This algorithm is based on the RV2COE algorithm presented in:
 * David A. Vallado, James Wertz - Fundamentals of Astrodynamics and Applications, 4th ed (Page 114)
 *
 */
template<typename Scalar>
std::tuple<Scalar, Scalar, Scalar, Scalar, Scalar, Scalar> rv2coe(
        Scalar mu, const Vector3<Scalar> &r, const Vector3<Scalar> &v, Scalar tol = 1e-8) {
    // h = r x v, h_mag = |h|
    Vector3<Scalar> h      = r.cross(v);
    Scalar          h_norm = h.norm();

    // n = k_hat x h, n_mag = |n|
    Vector3<Scalar> n      = Vector3<Scalar>::UnitZ().cross(h);
    Scalar          n_norm = n.norm();

    // e = ((v^2 - mu/r)*r - (r.v)*v)/mu, e_mag = |e|
    Vector3<Scalar> e      = (r * (v.squaredNorm() - mu / r.norm()) - v * r.dot(v)) / mu;
    Scalar          e_norm = e.norm();

    // p = h^2 / mu
    Scalar p = h.squaredNorm() / mu;

    // i = acos(h_z / |h|)
    Scalar inc = std::acos(h.z() / h_norm);

    // Special Cases
    bool circular   = e_norm < tol;
    bool equatorial = std::abs(inc) < tol;


    Scalar raan, argp, nu;

    if (equatorial && !circular) {
        raan = 0;
        argp = mod(std::atan2(e.y(), e.x()), 2 * odin::consts::PI<Scalar>);
        nu   = std::atan2(h.dot(e.cross(r)) / h_norm, e.dot(r));
    } else if (!equatorial && circular) {
        raan = mod(std::atan2(n.y(), n.x()), 2 * odin::consts::PI<Scalar>);
        argp = 0;
        nu   = std::atan2(r.dot(h.cross(n)) / h_norm, r.dot(n));
    } else if (equatorial && circular) {
        raan = 0;
        argp = 0;
        nu   = mod(std::atan2(r.y(), r.x()), 2 * odin::consts::PI<Scalar>);
    } else {
        Scalar a    = p / (1 - std::pow(e_norm, 2));
        Scalar mu_a = mu * a;
        if (a > 0) {
            Scalar e_se = r.dot(v) / std::sqrt(mu_a);
            Scalar e_ce = r.norm() * v.squaredNorm() / mu - 1;
            nu          = anomaly_eccentric_to_true(std::atan2(e_se, e_ce), e_norm);
        } else {
            Scalar e_sh = r.dot(v) / std::sqrt(-mu_a);
            Scalar e_ch = r.norm() * v.squaredNorm() / mu - 1;
            nu = anomaly_hyperbolic_to_true(std::log((e_ch + e_sh) / (e_ch - e_sh)) / 2, e_norm);
        }
        raan      = mod(std::atan2(n.y(), n.x()), 2 * odin::consts::PI<Scalar>);
        Scalar px = r.dot(n);
        Scalar py = r.dot(h.cross(n)) / h_norm;
        argp      = mod(std::atan2(py, px) - nu, 2 * odin::consts::PI<Scalar>);
    }
    nu = mod(nu + odin::consts::PI<Scalar>, 2 * odin::consts::PI<Scalar>) - odin::consts::PI<Scalar>;
    return std::make_tuple(p, e_norm, inc, raan, argp, nu);
}

//template<typename Scalar>
//constexpr auto rv2coe(Scalar                mu,
//        const std::vector<Vector3<Scalar>> &rr,
//        const std::vector<Vector3<Scalar>> &vv,
//        Scalar                              tol = 1e-8) {
//    std::vector<Scalar> p(rr.size());
//    std::vector<Scalar> e(rr.size());
//    std::vector<Scalar> inc(rr.size());
//    std::vector<Scalar> raan(rr.size());
//    std::vector<Scalar> argp(rr.size());
//    std::vector<Scalar> nu(rr.size());
//    for (size_t i = 0; i < rr.size(); ++i) {
//        std::tie(p[i], e[i], inc[i], raan[i], argp[i], nu[i]) = rv2coe(mu, rr[i], vv[i], tol);
//    }
//    return std::make_tuple(p, e, inc, raan, argp, nu);
//}

template<typename Scalar>
constexpr std::tuple<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> coe2pqr(
        Scalar p, Scalar e, Scalar nu, Scalar mu) {
    Scalar cos_nu    = std::cos(nu);
    Scalar sin_nu    = std::sin(nu);
    Scalar sqrt_mu_p = std::sqrt(mu / p);
    return std::make_tuple(
            Eigen::Vector3<Scalar>{p * cos_nu / (1 + e * cos_nu), p * sin_nu / (1 + e * cos_nu), 0},
            Eigen::Vector3<Scalar>{-sqrt_mu_p * sin_nu, sqrt_mu_p * (e + cos_nu), 0});
}

template<int N, typename Scalar>
constexpr Eigen::Matrix3<Scalar> rot(const Scalar &angle) {
    static_assert(N >= 0 && N <= 2, "N must be between 0 and 2 (inclusive and 0-indexed)");
    const Eigen::Vector3<Scalar> unit_vectors[] = {Eigen::Vector3<Scalar>::UnitX(),
            Eigen::Vector3<Scalar>::UnitY(),
            Eigen::Vector3<Scalar>::UnitZ()};
    return Eigen::AngleAxis<Scalar>(angle, unit_vectors[N]).toRotationMatrix();
}

template<typename Scalar>
constexpr std::tuple<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> coe2rv(
        Scalar mu,  // gravitational parameter
        Scalar p,   // semi-latus rectum
        Scalar e,   // eccentricity
        Scalar i,   // inclination
        Scalar raan,// right ascension of ascending node
        Scalar argp,// argument of periapsis
        Scalar nu   // true anomaly
) {
    // Need to implement the coe2pqr function, not provided in the original code
    auto [r_pqw, v_pqw] = coe2pqr<Scalar>(p, e, nu, mu);

    // create Eigen rotation matrices
    Eigen::Matrix3<Scalar> rot_ijk_pqw
            = rot<2, Scalar>(raan) * rot<0, Scalar>(i) * rot<2, Scalar>(argp);

    // return the rotated vectors
    return std::make_tuple(rot_ijk_pqw * r_pqw, rot_ijk_pqw * v_pqw);
}

//template<typename Scalar>
//constexpr std::tuple<Eigen::Matrix<Scalar, 3, Eigen::Dynamic>,
//        Eigen::Matrix<Scalar, 3, Eigen::Dynamic>>
//
//coe2rv(const Scalar                mu,  // gravitational parameter
//        const Scalar               p,   // semi-latus rectum
//        const Scalar               e,   // eccentricity
//        const Scalar               inc, // inclination
//        const Scalar               raan,// right ascension of ascending node
//        const Scalar               argp,// argument of periapsis
//        const std::vector<Scalar> &nu   // true anomaly
//) {
//    size_t                                   n = nu.size();
//    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> r_out(3, n);
//    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> v_out(3, n);
//
//#ifdef ODIN_USE_TBB
//    tbb::parallel_for(
//            tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range<size_t> &range) {
//                for (size_t i = range.begin(); i != range.end(); ++i) {
//                    auto temp    = coe2rv(mu, p, e, inc, raan, argp, nu[i]);
//                    r_out.col(i) = std::get<0>(temp);
//                    v_out.col(i) = std::get<1>(temp);
//                }
//            });
//#else
//    for (size_t i = 0; i < n; ++i) {
//        auto temp    = coe2rv(mu, p, e, inc, raan, argp, nu[i]);
//        r_out.col(i) = std::get<0>(temp);
//        v_out.col(i) = std::get<1>(temp);
//    }
//#endif
//    return std::make_tuple(r_out, v_out);
//}

//template<typename Scalar>
//constexpr std::tuple<Eigen::Matrix<Scalar, 3, Eigen::Dynamic>,
//        Eigen::Matrix<Scalar, 3, Eigen::Dynamic>>
//
//coe2rv(const Scalar                mu,  // gravitational parameter
//        const std::vector<Scalar> &p,   // semi-latus rectum
//        const std::vector<Scalar> &e,   // eccentricity
//        const std::vector<Scalar> &i,   // inclination
//        const std::vector<Scalar> &raan,// right ascension of ascending node
//        const std::vector<Scalar> &argp,// argument of periapsis
//        const std::vector<Scalar> &nu   // true anomaly
//) {
//    size_t                                   n = p.size();
//    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> r_out(3, n);
//    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> v_out(3, n);

//#ifdef ODIN_USE_TBB
//    tbb::parallel_for(tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range<size_t> &r) {
//        for (size_t idx = r.begin(); idx != r.end(); ++idx) {
//            auto [r_pqw, v_pqw]                = coe2pqr<Scalar>(p[idx], e[idx], nu[idx], mu);
//            Eigen::Matrix3<Scalar> rot_ijk_pqw = rot<2, Scalar>(raan[idx]) * rot<0, Scalar>(i[idx])
//                                               * rot<2, Scalar>(argp[idx]);
//            r_out.col(idx) = rot_ijk_pqw * r_pqw;
//            v_out.col(idx) = rot_ijk_pqw * v_pqw;
//        }
//    });
//#else
//    for (size_t idx = 0; idx < n; ++idx) {
//        auto [r_pqw, v_pqw] = coe2pqr<Scalar>(p[idx], e[idx], nu[idx], mu);
//        Eigen::Matrix3<Scalar> rot_ijk_pqw
//                = rot<2, Scalar>(raan[idx]) * rot<0, Scalar>(i[idx]) * rot<2, Scalar>(argp[idx]);
//        r_out.col(idx) = rot_ijk_pqw * r_pqw;
//        v_out.col(idx) = rot_ijk_pqw * v_pqw;
//    }
//#endif
//    return std::make_tuple(r_out, v_out);
//}

template<typename Scalar, typename ConversionFunc>
auto sample_anomalies(Scalar e, int n_points, ConversionFunc conversion) {
    Scalar step = 2 * odin::consts::PI<Scalar> / (n_points - 1);

    // Create a range, perform conversion
    // TODO: See if this is actually possible with LLVM 15/16. It's a bit of a hassle right now.
    //    auto anomalies = std::views::iota(0, n_points)
    //                   | std::views::transform([=](int i) { return i * step; })// Anomalies
    //                   | std::views::transform(conversion);                    // Perform conversion

    std::vector<Scalar> anomalies(n_points);
    for (int i = 0; i < n_points; i++) {
        Scalar anomaly = i * step;           // Anomalies
        anomalies[i]   = conversion(anomaly);// Perform conversion
    }
    return std::vector<Scalar>(anomalies.begin(), anomalies.end());
}

template<typename Scalar>
auto sample_true_from_mean_anomaly(Scalar e, int n_points) {
    // Define conversion from mean to true anomaly
    auto conversion = [=](Scalar m) {
        return anomaly_eccentric_to_true(anomaly_mean_to_eccentric(m, e), e);
    };
    return sample_anomalies(e, n_points, conversion);
}

template<typename Scalar>
auto sample_true_from_eccentric_anomaly(Scalar e, int n_points) {
    // Define conversion from eccentric to true anomaly
    auto conversion = [=](Scalar E) { return anomaly_eccentric_to_true(E, e); };
    return sample_anomalies(e, n_points, conversion);
}

}// namespace odin::domain::astrodynamics

#endif// ASTRODYNAMICS_HPP
