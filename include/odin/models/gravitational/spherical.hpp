#ifndef SPHERICAL_HPP
#define SPHERICAL_HPP

/*
 * @bibliography
 * - https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html
 * - https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html#Point-mass-potential
 * - https://phys.libretexts.org/Bookshelves/University_Physics/Book%3A_Physics_(Boundless)/5%3A_Uniform_Circular_Motion_and_Gravitation/5.5%3A_Newtons_Law_of_Universal_Gravitation
 * - https://en.wikipedia.org/wiki/Shell_theorem
 */
#include "gravitational_base.hpp"
#include <Eigen/Core>


template<typename Derived, typename S, size_t Dim = 3>
class SphericalBase {
public:
    using Scalar            = S;
    const static size_t dim = Dim;
    using Vector            = Eigen::Vector<Scalar, Dim>;

    explicit SphericalBase(const Scalar  mu,
                           const Scalar  radius   = 0,
                           const Vector &position = Vector::Zero())
        : position_(position), mu_(mu), radius_(radius) {}

    Scalar static_potential(const Vector &rel_position, const Scalar mu, const Scalar radius = 0) {
        return static_cast<Derived *>(this)->static_potential_impl(rel_position, mu, radius);
    }

    Vector static_acceleration(const Vector &rel_position, const Scalar mu, const Scalar radius = 0) {
        return static_cast<Derived *>(this)->static_acceleration_impl(rel_position, mu, radius);
    }

    [[nodiscard]] Scalar potential(const Vector &position) {
        return static_potential(position - position_, mu_, radius_);
    }

    [[nodiscard]] Vector acceleration(const Vector &position) {
        return static_acceleration(position - position_, mu_, radius_);
    }

    [[nodiscard]] SphericalBase thread_local_copy() const {
        SphericalBase copy = *this;
        return copy;
    }

protected:
    Vector position_;
    Scalar mu_;
    Scalar radius_;
};


template<typename S, size_t Dim = 3>
class PointMass : public SphericalBase<PointMass<S, Dim>, S, Dim> {
    using Base = SphericalBase<PointMass<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit PointMass(
            const Scalar  mu,
            const Vector &position = Vector::Zero())
        : Base(mu, 0, position) {}

    static Scalar static_potential_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar) {
        return -mu / rel_position.norm();
    }

    static Vector static_acceleration_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar) {
        return -mu * rel_position / pow(rel_position.norm(), 3);
    }
};


template<typename S, size_t Dim = 3>
class HollowSphere : public SphericalBase<HollowSphere<S, Dim>, S, Dim> {
    using Base = SphericalBase<HollowSphere<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit HollowSphere(
            const Scalar  mu,
            const Scalar  radius,
            const Vector &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        const Scalar r = rel_position.norm();
        return (r > radius)
                     ? PointMass<Scalar, Dim>::static_potential_impl(rel_position, mu, radius)
                     : -mu / radius;
    }
    static Vector static_acceleration_impl(
            const Vector &position,
            const Scalar  mu,
            const Scalar  radius) {
        const Scalar r = position.norm();
        return (r > radius)
                     ? PointMass<Scalar, Dim>::static_acceleration_impl(position, mu, radius)
                     : Vector::Zero();
    }
};

template<typename S, size_t Dim = 3>
class HomogeneousSphere : public SphericalBase<HomogeneousSphere<S, Dim>, S, Dim> {
    using Base = SphericalBase<HomogeneousSphere<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit HomogeneousSphere(
            const Scalar  mu,
            const Scalar  radius,
            const Vector &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        const Scalar r = rel_position.norm();
        return (r > radius)
                     ? PointMass<Scalar, Dim>::static_potential_impl(rel_position, mu, radius)
                     : -mu * (3 * radius * radius - r * r) / (2 * pow(radius, 3));
    }
    static Vector static_acceleration_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        const Scalar r = rel_position.norm();
        return (r > radius)
                     ? PointMass<Scalar, Dim>::static_acceleration_impl(rel_position, mu, radius)
                     : -mu * rel_position / pow(radius, 3);
    }
};

template<typename S, size_t Dim = 3>
class PlummerSphere : public SphericalBase<PlummerSphere<S, Dim>, S, Dim> {
    using Base = SphericalBase<PlummerSphere<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit PlummerSphere(
            const Scalar  mu,
            const Scalar  radius,
            const Vector &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        return -mu / std::sqrt(rel_position.squaredNorm() + radius * radius);
    }

    static Vector static_acceleration_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        return -mu * rel_position.normalized() / (rel_position.squaredNorm() + radius * radius);
    }
};

template<typename S, size_t Dim = 3>
class [[maybe_unused]] IsochroneSphere : public SphericalBase<IsochroneSphere<S, Dim>, S, Dim> {
    using Base = SphericalBase<IsochroneSphere<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit IsochroneSphere(
            const Scalar  mu,
            const Scalar  radius,
            const Vector &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        const Scalar r = rel_position.norm();
        return -mu / (radius + std::sqrt(r * r + radius * radius));
    }

    static Vector static_acceleration_impl(
            const Vector &rel_position,
            const Scalar  mu,
            const Scalar  radius) {
        const Scalar r = rel_position.norm();
        return -mu * rel_position / pow(radius + std::sqrt(r * r + radius * radius), 3);
    }
};


template<typename Scalar, int Degree, int Order>
struct CoefficientsHolder {
    using Clm = std::array<std::array<Scalar, Order + 1>, Degree + 1>;
    using Slm = std::array<std::array<Scalar, Order + 1>, Degree + 1>;
};

template<typename Scalar>
struct CoefficientsHolder<Scalar, -1, -1> {
    using Clm = std::vector<std::vector<Scalar>>;
    using Slm = std::vector<std::vector<Scalar>>;
};

template<typename Scalar, int Degree, int Order>
struct AssociatedLegendrePolynomialsHolder {
    using Plm = std::array<std::array<Scalar, Order + 1>, Degree + 1>;
    static Plm init() {
        return Plm{};
    }
};

template<typename Scalar>
struct AssociatedLegendrePolynomialsHolder<Scalar, -1, -1> {
    using Plm = std::vector<std::vector<Scalar>>;
    static Plm init(int degree, int order) {
        return Plm(degree + 1, std::vector<Scalar>(order + 1));
    }
};

// Standard factorial function
//unsigned long long factorial(int n) {
//    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
//}

template<typename Scalar>
Scalar factorial(int n) {
    Scalar result = 1;

    for (int i = 1; i <= n; ++i) {
        result *= i;
    }

    return result;
}

template<typename Scalar, typename Plm, int Degree = -1, int Order = -1>
Plm normalize_legendre(const Plm &P_lm) {
    Plm Pbar_lm;

    // Make sure to resize for dynamic containers
    if constexpr (Degree < 0 && Order < 0) {
        Pbar_lm.resize(P_lm.size(), std::vector<Scalar>(P_lm[0].size()));
    }

    for (int l = 0; l < P_lm.size(); l++) {
        for (int m = 0; m < P_lm[l].size(); m++) {
            Scalar factor = sqrt((2 - (m == 0 ? 1 : 0)) * (2 * l + 1) * factorial<Scalar>(l - m) / factorial<Scalar>(l + m));
            Pbar_lm[l][m] = factor * P_lm[l][m];
        }
    }

    return Pbar_lm;
}

template<typename Scalar, typename Plm, int Degree = -1, int Order = -1>
Plm denormalize_legendre(const Plm &Pbar_lm) {
    Plm P_lm;

    // Make sure to resize for dynamic containers
    if constexpr (Degree < 0 && Order < 0) {
        P_lm.resize(Pbar_lm.size(), std::vector<Scalar>(Pbar_lm[0].size()));
    }

    for (int l = 0; l < Pbar_lm.size(); l++) {
        for (int m = 0; m < Pbar_lm[l].size(); m++) {
            Scalar factor = sqrt(1 / ((2 - (m == 0 ? 1 : 0)) * (2 * l + 1)) * factorial<Scalar>(l + m) / factorial<Scalar>(l - m));
            P_lm[l][m]    = factor * Pbar_lm[l][m];
        }
    }

    return P_lm;
}


//// Base class for common functions
//template<typename Derived, typename Scalar>
//struct CoefficientTransformationsBase {
//    using Clm = typename Derived::Clm;
//    using Slm = typename Derived::Slm;
//
//
//};

// General template for any Degree and Order
template<typename Scalar, int Degree, int Order>
struct CoefficientTransformations {
    using Clm = typename CoefficientsHolder<Scalar, Degree, Order>::Clm;
    using Slm = typename CoefficientsHolder<Scalar, Degree, Order>::Slm;

    std::pair<Clm, Slm> normalized(const Clm &inputClm, const Slm &inputSlm) {
        Clm normalizedClm = inputClm;
        Slm normalizedSlm = inputSlm;
        normalize_in_place(normalizedClm, normalizedSlm);
        return {normalizedClm, normalizedSlm};
    }

    std::pair<Clm, Slm> denormalized(const Clm &inputClm, const Slm &inputSlm) {
        Clm denormalizedClm = inputClm;
        Slm denormalizedSlm = inputSlm;
        denormalize_in_place(denormalizedClm, denormalizedSlm);
        return {denormalizedClm, denormalizedSlm};
    }

    static void normalize_in_place(Clm &inputClm, Slm &inputSlm) {
        for (int l = 0; l <= Degree; ++l) {
            for (int m = 0; m <= std::min(l, Order); ++m) {
                Scalar factor = sqrt(1.0 / ((2 - (m == 0 ? 1 : 0)) * (2 * l + 1)) * factorial<Scalar>(l + m) / factorial<Scalar>(l - m));
                inputClm[l][m] *= factor;
                inputSlm[l][m] *= factor;
            }
        }
    }

    static void denormalize_in_place(Clm &inputClm, Slm &inputSlm) {
        for (int l = 0; l <= Degree; ++l) {
            for (int m = 0; m <= std::min(l, Order); ++m) {
                Scalar factor = sqrt((2 - (m == 0 ? 1 : 0)) * (2 * l + 1) * factorial<Scalar>(l - m) / factorial<Scalar>(l + m));
                inputClm[l][m] *= factor;
                inputSlm[l][m] *= factor;
            }
        }
    }
};

// Partial specialization for dynamic Degree and Order
template<typename Scalar>
struct CoefficientTransformations<Scalar, -1, -1> {
    // Implementation for dynamic Degree and Order
    using Clm = typename CoefficientsHolder<Scalar, -1, -1>::Clm;
    using Slm = typename CoefficientsHolder<Scalar, -1, -1>::Slm;

    std::pair<Clm, Slm> normalized(const Clm &inputClm, const Slm &inputSlm) {
        Clm normalizedClm = inputClm;
        Slm normalizedSlm = inputSlm;
        normalize_in_place(normalizedClm, normalizedSlm);
        return {normalizedClm, normalizedSlm};
    }

    std::pair<Clm, Slm> denormalized(const Clm &inputClm, const Slm &inputSlm) {
        Clm denormalizedClm = inputClm;
        Slm denormalizedSlm = inputSlm;
        denormalize_in_place(denormalizedClm, denormalizedSlm);
        return {denormalizedClm, denormalizedSlm};
    }

    static void normalize_in_place(Clm &inputClm, Slm &inputSlm) {
        for (int l = 0; l < inputClm.size(); ++l) {
            for (int m = 0; m < inputClm[l].size(); ++m) {
                Scalar factor = sqrt(1.0 / ((2 - (m == 0 ? 1 : 0)) * (2 * l + 1)) * factorial<Scalar>(l + m) / factorial<Scalar>(l - m));
                inputClm[l][m] *= factor;
                inputSlm[l][m] *= factor;
            }
        }
    }

    static void denormalize_in_place(Clm &inputClm, Slm &inputSlm) {
        for (int l = 0; l < inputClm.size(); ++l) {
            for (int m = 0; m < inputClm[l].size(); ++m) {
                Scalar factor = sqrt((2 - (m == 0 ? 1 : 0)) * (2 * l + 1) * factorial<Scalar>(l - m) / factorial<Scalar>(l + m));
                inputClm[l][m] *= factor;
                inputSlm[l][m] *= factor;
            }
        }
    }
};


template<typename Scalar, int Degree, int Order>
struct GeopotentialCoefficientTransformations {
    // Implementation for dynamic Degree and Order
    using Clm = typename CoefficientsHolder<Scalar, Degree, Order>::Clm;
    using Slm = typename CoefficientsHolder<Scalar, Degree, Order>::Slm;

    static std::pair<Clm, Slm> add_dimensions(Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
        Clm C_lm = inputC_lm;
        Slm S_lm = inputS_lm;
        add_dimensions_in_place(mu, R, C_lm, S_lm);
        return {C_lm, S_lm};
    }

    std::pair<Clm, Slm> remove_dimensions(Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
        Clm C_lm = inputC_lm;
        Slm S_lm = inputS_lm;
        remove_dimensions_in_place(mu, R, C_lm, S_lm);
        return {C_lm, S_lm};
    }

    static void add_dimensions_in_place(Scalar mu, Scalar R, Clm &C_lm, Slm &S_lm) {
        for (int l = 0; l <= Degree; l++) {
            Scalar scale_factor = -mu * std::pow(R, l);
            for (int m = 0; m <= std::min(l, Order); m++) {
                C_lm[l][m] *= scale_factor;
                S_lm[l][m] *= scale_factor;
            }
        }
    }

    static void remove_dimensions_in_place(Scalar mu, Scalar R, Clm &C_lm, Slm &S_lm) {
        for (int l = 0; l <= Degree; l++) {
            Scalar scale_factor = -mu * std::pow(R, -l);
            for (int m = 0; m <= std::min(l, Order); m++) {
                C_lm[l][m] *= scale_factor;
                S_lm[l][m] *= scale_factor;
            }
        }
    }
};

// Partial specialization for dynamic Degree and Order
template<typename Scalar>
struct GeopotentialCoefficientTransformations<Scalar, -1, -1> {
    using Clm = typename CoefficientsHolder<Scalar, -1, -1>::Clm;
    using Slm = typename CoefficientsHolder<Scalar, -1, -1>::Slm;

    std::pair<Clm, Slm> add_dimensions(Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
        Clm C_lm = inputC_lm;
        Slm S_lm = inputS_lm;
        add_dimensions_in_place(mu, R, C_lm, S_lm);
        return {C_lm, S_lm};
    }

    std::pair<Clm, Slm> remove_dimensions(Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
        Clm C_lm = inputC_lm;
        Slm S_lm = inputS_lm;
        remove_dimensions_in_place(mu, R, C_lm, S_lm);
        return {C_lm, S_lm};
    }

    static void add_dimensions_in_place(Scalar mu, Scalar R, Clm &C_lm, Slm &S_lm) {
        for (int l = 0; l < C_lm.size(); l++) {
            Scalar scale_factor = -mu * std::pow(R, l);
            for (int m = 0; m < C_lm[l].size(); m++) {
                C_lm[l][m] *= scale_factor;
                S_lm[l][m] *= scale_factor;
            }
        }
    }
    static void remove_dimensions_in_place(Scalar mu, Scalar R, Clm &C_lm, Slm &S_lm) {
        for (int l = 0; l < C_lm.size(); l++) {
            Scalar scale_factor = 1 / (-mu * std::pow(R, l));
            for (int m = 0; m < C_lm[l].size(); m++) {
                C_lm[l][m] *= scale_factor;
                S_lm[l][m] *= scale_factor;
            }
        }
    }
};


//// Compile-time
//CoefficientsHolder<double, 2, 2>::Clm Clm = {...};
//CoefficientsHolder<double, 2, 2>::Slm Slm = {...};
//auto [normalizedClm, normalizedSlm] = ForwardNormalize<double, 2, 2>::normalize(Clm, Slm);
//
//// Runtime
//CoefficientsHolder<double, -1, -1>::Clm Clm = {...};
//CoefficientsHolder<double, -1, -1>::Slm Slm = {...};
//auto [normalizedClm, normalizedSlm] = ForwardNormalize<double, -1, -1>::normalize(Clm, Slm);


//// Compile-time
//auto [denormalizedClm, denormalizedSlm] = BackwardNormalize<double, 2, 2>::denormalize(normalizedClm, normalizedSlm);
//
//// Runtime
//auto [denormalizedClm, denormalizedSlm] = BackwardNormalize<double, -1, -1>::denormalize(normalizedClm, normalizedSlm);

template<typename S, int Degree, int Order = Degree, int Dim = 3>
class [[maybe_unused]] SphericalHarmonics : public SphericalBase<SphericalHarmonics<S, Degree, Order, Dim>, S, Dim> {
    using Base = SphericalBase<SphericalHarmonics<S, Degree, Order, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = S;
    using Clm    = typename CoefficientsHolder<Scalar, Degree, Order>::Clm;
    using Slm    = typename CoefficientsHolder<Scalar, Degree, Order>::Slm;
    using Plm    = typename AssociatedLegendrePolynomialsHolder<Scalar, Degree, Order>::Plm;
    Clm Clm_;
    Slm Slm_;
    Plm Plm_;

    explicit SphericalHarmonics(
            const Scalar  mu,
            const Scalar  radius,
            const Clm     Cnm      = CoefficientsHolder<Scalar, Degree, Order>::Clm(),
            const Slm     Snm      = CoefficientsHolder<Scalar, Degree, Order>::Slm(),
            const Vector &position = Vector::Zero())
        : Base(mu, radius, position), Clm_(Cnm), Slm_(Snm) {
        if constexpr (Degree == -1 && Order == -1) {
            int degree = Clm_.size() - 1;
            int order  = Clm_.front().size() - 1;
            Plm_       = AssociatedLegendrePolynomialsHolder<Scalar, Degree, Order>::init(degree, order);
        } else {
            Plm_ = AssociatedLegendrePolynomialsHolder<Scalar, Degree, Order>::init();
        }
    }

    void set_coefficients(const Clm &Cnm, const Slm &Snm) {
        Clm_ = Cnm;
        Slm_ = Snm;
    }

    Scalar static_potential_impl(const Vector &rel_position, const Scalar mu, const Scalar radius) {
        Scalar r = rel_position.norm();
        if (r < radius) {
            return HomogeneousSphere<Scalar, Dim>::static_potential_impl(rel_position, mu, radius);
        }

        Scalar phi     = std::acos(rel_position.z() / r);
        Scalar lambda  = std::atan2(rel_position.y(), rel_position.x());
        Scalar sin_phi = std::sin(phi);

        Scalar V               = 0;
        Scalar a_r_ratio_pow_l = 1;// to keep track of (a_e / r) ^ l

        for (size_t l = 0; l < Clm_.size(); ++l) {
            for (size_t m = 0; m <= l && m < Clm_[l].size(); ++m) {
                Scalar cos_m_lambda = std::cos(m * lambda);
                Scalar sin_m_lambda = std::sin(m * lambda);
                Scalar P_lm         = std::assoc_legendre(l, m, sin_phi);
                Scalar term         = P_lm * (Clm_[l][m] * cos_m_lambda + Slm_[l][m] * sin_m_lambda);
                V += a_r_ratio_pow_l * term;
                //                ODIN_LOG_INFO << "l: " << l << " m: " << m << " P_lm: " << P_lm << " term: " << term << " V: " << V << std::endl;
            }
            a_r_ratio_pow_l *= radius / r;
        }
        V *= mu / r;

        return V;
    }

    constexpr Vector spherical_to_cartesian(const Vector &spherical, Scalar lambda, Scalar phi) {
        Scalar                 sin_lambda = std::sin(lambda);
        Scalar                 cos_lambda = std::cos(lambda);
        Scalar                 cos_phi    = std::cos(phi);
        Scalar                 sin_phi    = std::sin(phi);
        Eigen::Matrix3<Scalar> transform{
                {cos_phi * cos_lambda, -sin_phi * cos_lambda, -sin_lambda},
                {cos_phi * sin_lambda, -sin_phi * sin_lambda,  cos_lambda},
                {             sin_phi,               cos_phi,           0}
        };
        return transform * spherical;
    }


    Vector static_acceleration_impl(const Vector &rel_position, const Scalar mu, const Scalar radius) {
        Scalar r = rel_position.norm();
        if (r < radius) {
            return HomogeneousSphere<Scalar, Dim>::static_acceleration_impl(rel_position, mu, radius);
        }

        Scalar phi    = std::acos(rel_position.z() / r);
        Scalar lambda = std::atan2(rel_position.y(), rel_position.x());

        Scalar sin_phi         = std::sin(phi);
        Scalar cos_phi         = std::cos(phi);
        Scalar r2_inv          = -mu / (r * r);
        Scalar a_r_ratio_pow_l = 1;// to keep track of (a_e / r) ^ l

        Vector acceleration = Vector::Zero();

        //        //
        //        Plm_[0][0] = 1;
        //        //
        //        for (size_t l = 1; l < Clm_.size(); ++l) {
        //            Scalar P_ll = std::assoc_legendre(l, l, sin_phi);
        //            Scalar P_ll_1;
        //            if (sin_phi == 1 || sin_phi == -1) {
        //                P_ll_1 = 0;// or some other appropriate value
        //            } else {
        //                P_ll_1 = ((2 * l - 1) * sin_phi * P_ll - (l - 1) * Plm_[l - 1][l - 1]) / sqrt(1 - sin_phi * sin_phi);
        //            }
        //            Plm_[l][l]     = P_ll;
        //            Plm_[l][l - 1] = P_ll_1;
        //        }
        //        //
        //        //        // now start from l=2 and go down the table
        //        for (size_t l = 2; l < Clm_.size(); ++l) {
        //            for (size_t m = 0; m < l - 1; ++m) {
        //                Plm_[l][m] = ((2 * l - 1) * sin_phi * Plm_[l - 1][m] - (l + m - 1) * Plm_[l - 2][m]) / (l - m);
        //            }
        //        }

        // precompute all P
        for (size_t l = 0; l < Clm_.size(); ++l) {
            for (size_t m = 0; m <= l; ++m) {
                Plm_[l][m] = std::assoc_legendre(l, m, sin_phi);
            }
        }

        for (size_t l = 0; l < Clm_.size(); ++l) {
            for (size_t m = 0; m <= l; ++m) {
                Scalar cos_m_lambda = std::cos(m * lambda);
                Scalar sin_m_lambda = std::sin(m * lambda);
                Scalar common_term  = (Clm_[l][m] * cos_m_lambda + Slm_[l][m] * sin_m_lambda);
                acceleration[0] -= ((l + 1) * Plm_[l][m] * common_term) * a_r_ratio_pow_l;// radial component

                if (l == 0 && m == 0) { continue; }                                                                     // only radial component for l=0, m=0
                Scalar x = sin_phi;  // We make the change of variables x = sin(phi)
                Scalar dx_dphi = cos_phi;  // Derivative of sin(phi) w.r.t. phi

                // Compute the derivative of the associated Legendre function w.r.t. x
                Scalar dPlm_dx = (1- x * x) * ((l - m + 1) * x * Plm_[l][m] - (l + 1) * Plm_[l + 1][m]);

                // Then, by the chain rule, compute the derivative of P_{l, m}(sin(phi)) w.r.t. phi
                Scalar dPlm_dphi = dPlm_dx * dx_dphi;
                acceleration[1] += (dPlm_dphi * common_term) * a_r_ratio_pow_l;

                if (l == 1 && m == 0) { continue; }
                acceleration[2] += m * (Plm_[l][m] / cos_phi * (-Clm_[l][m] * sin_m_lambda + Slm_[l][m] * cos_m_lambda)) * a_r_ratio_pow_l;
            }
            a_r_ratio_pow_l *= radius / r;
        }
        acceleration *= r2_inv;
        return spherical_to_cartesian(acceleration, lambda, phi);
    }
};

// Verify that it adheres to the GravitationalModel concept
static_assert(is_gravitational_model_v<PointMass<double>>, "HollowSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<HollowSphere<double>>, "HollowSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<HomogeneousSphere<double>>, "HomogeneousSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<PlummerSphere<double>>, "PlummerSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<IsochroneSphere<double>>, "IsochroneSphere does not comply with GravitationalModel.");

#endif// SPHERICAL_HPP
