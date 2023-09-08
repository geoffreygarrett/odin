#ifndef SPHERICAL_HPP
#define SPHERICAL_HPP

/*
 * @bibliography
 * - https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html
 * - https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html#Point-mass-potential
 * - https://phys.libretexts.org/Bookshelves/University_Physics/Book%3A_Physics_(Boundless)/5%3A_Uniform_Circular_Motion_and_Gravitation/5.5%3A_Newtons_Law_of_Universal_Gravitation
 * - https://en.wikipedia.org/wiki/Shell_theorem
 */
#include <Eigen/Core>

#include "gravitational_base.hpp"
#include "gravitational_concept.hpp"

#ifdef ODIN_AUTODIFF
#define ODIN_CONST
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
using namespace autodiff;
#else
#define ODIN_CONST const
#endif

namespace mlib {


#ifdef __APPLE__// if this is apple, use custom implementation of assoc_legendre
#include <cmath>

template<typename Scalar>
Scalar factorial(int n) {
    if (n == 0) return 1.0;
    Scalar result = 1.0;
    for (int i = 1; i <= n; ++i) { result *= i; }
    return result;
}

template<typename Scalar>
Scalar assoc_legendre(int l, int m, double x) {
    Scalar sum = 0.0;
    for (int k = 0; k <= (l - m) / 2; ++k) {
        Scalar term = std::pow(-1.0, k) * factorial<Scalar>(2 * l - 2 * k)
                    / (std::pow(2.0, l) * factorial<Scalar>(k) * factorial<Scalar>(l - k)
                            * factorial<Scalar>(l - m - 2 * k))
                    * std::pow(1.0 - x * x, (l - m - 2 * k) / 2.0);
        sum += term;
    }
    return std::sqrt((2.0 * l + 1.0) * factorial<Scalar>(l - m) / (2.0 * factorial<Scalar>(l + m)))
         * sum;
}
#else

#include <cmath>

using assoc_legendre = std::assoc_legendre;

#endif
}// namespace mlib


// macro to check for C++20
#if __cplusplus > 201703L
// Two argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_2(TYPE1, TYPE2) requires std::is_same_v<TYPE1, TYPE2>
// Three argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_3(TYPE1, TYPE2, TYPE3) \
    requires(std::is_same_v<TYPE1, TYPE2> || std::is_same_v<TYPE1, TYPE3>)
#else
// Two argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_2(TYPE1, TYPE2) \
    std::enable_if_t<std::is_same_v<TYPE1, TYPE2>, TYPE2>
// Three argument version
#define ODIN_ENABLE_IF_TYPE_MATCHES_3(TYPE1, TYPE2, TYPE3)                                       \
    std::enable_if_t<std::disjunction_v<std::is_same<TYPE1, TYPE2>, std::is_same<TYPE1, TYPE3>>, \
            TYPE1>
#endif

// Using helper macros to differentiate between the number of arguments
#define GET_MACRO(_1, _2, _3, NAME, ...) NAME
#define ODIN_ENABLE_IF_TYPE_MATCHES(...)                                                 \
    GET_MACRO(__VA_ARGS__, ODIN_ENABLE_IF_TYPE_MATCHES_3, ODIN_ENABLE_IF_TYPE_MATCHES_2) \
    (__VA_ARGS__)


template<typename Derived, typename Scalar, size_t Dim = 3, typename... Params>
class model_base : public odin::crtp_base<model_base<Derived, Scalar, Dim, Params...>> {
public:
    using params_type = std::tuple<Params...>;
    using scalar_type = Scalar;
    using vector_type = Eigen::Matrix<scalar_type, Dim, 1>;

    explicit model_base(Params... params)
        : m_params(params...) {}

    // With eval
    template<typename Func, typename O, typename EvalType, typename... WrtArgs, typename... AtArgs>
    ODIN_ENABLE_IF_TYPE_MATCHES(scalar_type, autodiff::dual, autodiff::real)
    auto calc_partial(Func         lambda,
            EvalType              &eval,
            std::tuple<WrtArgs...> wrt_args,
            std::tuple<AtArgs...>  at_args) {
        if constexpr (std::is_same_v<O, scalar_type>) {
            return gradient(lambda,
                    wrt(std::get<WrtArgs>(wrt_args)...),
                    at(std::get<AtArgs>(at_args)...),
                    eval);
        } else if constexpr (std::is_same_v<O, vector_type>) {
            return jacobian(lambda,
                    wrt(std::get<WrtArgs>(wrt_args)...),
                    at(std::get<AtArgs>(at_args)...),
                    eval);
        } else {
            static_assert(sizeof(O) == -1, "Invalid return type for partial deducing template");
        }
    }

    // Without eval
    template<typename Func, typename O, typename... WrtArgs, typename... AtArgs>
    ODIN_ENABLE_IF_TYPE_MATCHES(scalar_type, autodiff::dual, autodiff::real)
    auto calc_partial(Func lambda, std::tuple<WrtArgs...> wrt_args, std::tuple<AtArgs...> at_args) {
        if constexpr (std::is_same_v<O, scalar_type>) {
            scalar_type eval;
            return gradient(lambda,
                    wrt(std::get<WrtArgs>(wrt_args)...),
                    at(std::get<AtArgs>(at_args)...),
                    eval);
        } else if constexpr (std::is_same_v<O, vector_type>) {
            vector_type eval;
            return jacobian(lambda,
                    wrt(std::get<WrtArgs>(wrt_args)...),
                    at(std::get<AtArgs>(at_args)...),
                    eval);
        } else {
            static_assert(sizeof(O) == -1, "Invalid return type for partial deducing template");
        }
    }

    // Function to calculate partial w.r.t. a subset of the parameters
    template<typename Func, typename O, std::size_t... I>
    auto calc_partial_with_params(Func lambda, std::index_sequence<I...>) {
        return calc_partial<Func, O>(
                [&](auto &&...args) {
                    return lambda(std::get<I>(m_params)..., std::forward<decltype(args)>(args)...);
                },
                std::get<I>(m_params)...);
    }

    // Function to calculate partial w.r.t. a single parameter at index I
    template<typename Func, typename O, std::size_t I>
    auto calc_partial_with_param(Func lambda) {
        return calc_partial<Func, O>(
                [&](auto &&...args) {
                    return lambda(std::get<I>(m_params), std::forward<decltype(args)>(args)...);
                },
                std::get<I>(m_params));
    }

protected:
    params_type m_params;
};


// Using named constants for better code readability
constexpr std::size_t CENTRAL_POSITION  = 0;
constexpr std::size_t GRAVITY_PARAMETER = 1;

template<typename Derived, typename Scalar, std::size_t Dim = 3>
class gravity_base
    : public model_base<Derived, Scalar, Dim, Scalar, Eigen::Matrix<Scalar, Dim, 1>> {
public:
    using vector_type = Eigen::Matrix<Scalar, Dim, 1>;
    using model_base<Derived, Scalar, Dim, Scalar, vector_type>::model_base;

    // Thread-local copy of the derived class
    Derived thread_local_copy() {
        return this->as_derived().thread_local_copy_impl();
    }

    // Potential energy calculation
    Scalar potential(ODIN_CONST vector_type &position) const {
        return this->as_derived().potential_impl(position);
    }

    // Acceleration calculation
    vector_type acceleration(ODIN_CONST vector_type &position) const {
        return this->as_derived().acceleration_impl(position);
    }
};

template<typename Scalar, std::size_t Dim = 3>
class point_mass_gravity : public gravity_base<point_mass_gravity<Scalar, Dim>, Scalar, Dim> {
public:
    using scalar_type = Scalar;
    using vector_type = Eigen::Matrix<Scalar, Dim, 1>;
    using gravity_base<point_mass_gravity<Scalar, Dim>, Scalar, Dim>::gravity_base;

    template<typename T>
    auto thread_local_copy_impl() {
        return point_mass_gravity<T, Dim>(std::get<CENTRAL_POSITION>(this->params_),
                std::get<GRAVITY_PARAMETER>(this->params_));
    }

    scalar_type potential_impl(ODIN_CONST vector_type &position) const {
        auto r  = position - std::get<CENTRAL_POSITION>(this->params_);
        auto GM = std::get<GRAVITY_PARAMETER>(this->params_);
        return -GM / r.norm();
    }

    vector_type acceleration_impl(ODIN_CONST vector_type &position) const {
        auto r           = position - std::get<CENTRAL_POSITION>(this->params_);
        auto GM          = std::get<GRAVITY_PARAMETER>(this->params_);
        auto r_norm_cube = r.norm() * r.norm() * r.norm();
        return -GM * r / r_norm_cube;
    }
};

template<typename Derived, typename S, size_t Dim = 3, typename... Params>
class SphericalBase {
public:
    using Scalar            = S;
    const static size_t dim = Dim;
    using Vector            = Eigen::Vector<Scalar, Dim>;

    explicit SphericalBase(ODIN_CONST Scalar mu,
            ODIN_CONST Scalar                radius   = 0,
            ODIN_CONST Vector               &position = Vector::Zero())
        : position_(position),
          mu_(mu),
          radius_(radius) {}

    Scalar static_potential(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius = 0) {
        return static_cast<Derived *>(this)->static_potential_impl(rel_position, mu, radius);
    }

    Vector static_acceleration(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius = 0) {
        return static_cast<Derived *>(this)->static_acceleration_impl(rel_position, mu, radius);
    }

#ifdef ODIN_AUTODIFF

    [[nodiscard]] Scalar potential(ODIN_CONST Vector &position) {
        Vector temp = position - position_;
        return static_potential(temp, mu_, radius_);
    }

    [[nodiscard]] Vector acceleration(ODIN_CONST Vector &position) {
        Vector temp = position - position_;
        return static_acceleration(temp, mu_, radius_);
    }

    [[nodiscard]] auto potential(ODIN_CONST Vector &&position) {
        Vector position_copy = std::move(position);// Create a non-const local copy
        Vector rel_position  = (position_copy - position_).eval();
        return static_potential(rel_position, mu_, radius_);
    }

    [[nodiscard]] auto acceleration(ODIN_CONST Vector &&position) {
        Vector position_copy = std::move(position);// Create a non-const local copy
        Vector rel_position  = (position_copy - position_).eval();
        return static_acceleration(rel_position, mu_, radius_);
    }

    template<typename Func, typename O, typename T = Scalar>
    ODIN_ENABLE_IF_TYPE_MATCHES(T, autodiff::dual, autodiff::real)
    auto calc_partial(Func lambda, ODIN_CONST Vector &rel_position) {
        if constexpr (std::is_same_v<O, Scalar>) {
            T eval;
            return gradient(lambda, wrt(rel_position), at(rel_position), eval);
        } else if constexpr (std::is_same_v<O, Vector>) {
            Vector eval;
            return jacobian(lambda, wrt(rel_position), at(rel_position), eval);
        } else {
            static_assert(sizeof(O) == -1, "Invalid return type for partial deducing template");
        }
    }

    // For lvalue references
    template<typename T = Scalar>
    ODIN_ENABLE_IF_TYPE_MATCHES(T, autodiff::dual, autodiff::real)
    auto partial_potential_wrt_position(ODIN_CONST Vector &position) {
        Vector rel_position = (position - position_).eval();
        auto   lambda       = [this, rel_position](ODIN_CONST auto &arg) {
            return this->static_potential(arg, this->mu_, this->radius_);
        };
        return calc_partial<decltype(lambda), Scalar, T>(lambda, rel_position);
    }

    // For rvalue references
    template<typename T = Scalar>
    ODIN_ENABLE_IF_TYPE_MATCHES(T, autodiff::dual, autodiff::real)
    auto partial_potential_wrt_position(ODIN_CONST Vector &&position) {
        Vector position_copy = std::move(position);
        return partial_potential_wrt_position<T>(position_copy);
    }

    // Partial acceleration wrt position for lvalue references
    template<typename T = Scalar>
    ODIN_ENABLE_IF_TYPE_MATCHES(T, autodiff::dual, autodiff::real)
    auto partial_acceleration_wrt_position(ODIN_CONST Vector &position) {
        Vector rel_position = (position - position_).eval();
        auto   lambda       = [this, &rel_position, mu_ = this->mu_, radius_ = this->radius_](
                              ODIN_CONST auto &rel_position_) {
            return this->static_acceleration(rel_position_, mu_, radius_);
        };
        return calc_partial<decltype(lambda), Vector, T>(lambda, rel_position);
    }

    // Partial acceleration wrt position for current instance, for rvalue references
    template<typename T = Scalar>
    ODIN_ENABLE_IF_TYPE_MATCHES(T, autodiff::dual, autodiff::real)
    auto partial_acceleration_wrt_position(ODIN_CONST Vector &&position) {
        Vector position_copy = std::move(position);             // Create a non-const local copy
        return partial_acceleration_wrt_position(position_copy);// Delegate to lvalue overload
    }

#else

    [[nodiscard]] Scalar potential(ODIN_CONST Vector &position) {
        return static_potential(position - position_, mu_, radius_);
    }

    [[nodiscard]] Vector acceleration(ODIN_CONST Vector &position) {
        return static_acceleration(position - position_, mu_, radius_);
    }

#endif

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

#ifndef ODIN_AUTODIFF
    explicit PointMass(ODIN_CONST Scalar mu, ODIN_CONST Vector &position = Vector::Zero())
#else
    explicit PointMass(ODIN_CONST Scalar mu, ODIN_CONST Vector &&position = Vector::Zero())
#endif
        : Base(mu, 0, position) {
    }

    static Scalar static_potential_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar) {
        return -mu / rel_position.norm();
    }

    static Vector static_acceleration_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar) {
        return -mu * rel_position / pow(rel_position.norm(), 3);
    }
};


template<typename S, size_t Dim = 3>
class HollowSphere : public SphericalBase<HollowSphere<S, Dim>, S, Dim> {
    using Base = SphericalBase<HollowSphere<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit HollowSphere(ODIN_CONST Scalar mu,
            ODIN_CONST Scalar               radius,
            ODIN_CONST Vector              &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        ODIN_CONST Scalar r = rel_position.norm();
        return (r > radius)
                     ? PointMass<Scalar, Dim>::static_potential_impl(rel_position, mu, radius)
                     : -mu / radius;
    }
    static Vector static_acceleration_impl(
            ODIN_CONST Vector &position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        ODIN_CONST Scalar r = position.norm();
        return (r > radius) ? PointMass<Scalar, Dim>::static_acceleration_impl(position, mu, radius)
                            : Vector::Zero();
    }
};

template<typename S, size_t Dim = 3>
class HomogeneousSphere : public SphericalBase<HomogeneousSphere<S, Dim>, S, Dim> {
    using Base = SphericalBase<HomogeneousSphere<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit HomogeneousSphere(ODIN_CONST Scalar mu,
            ODIN_CONST Scalar                    radius,
            ODIN_CONST Vector                   &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        ODIN_CONST Scalar r = rel_position.norm();
        return (r > radius)
                     ? PointMass<Scalar, Dim>::static_potential_impl(rel_position, mu, radius)
                     : -mu * (3 * radius * radius - r * r) / (2 * pow(radius, 3));
    }
    static Vector static_acceleration_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        ODIN_CONST Scalar r = rel_position.norm();
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

    explicit PlummerSphere(ODIN_CONST Scalar mu,
            ODIN_CONST Scalar                radius,
            ODIN_CONST Vector               &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        return -mu / std::sqrt(rel_position.squaredNorm() + radius * radius);
    }

    static Vector static_acceleration_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        return -mu * rel_position.normalized() / (rel_position.squaredNorm() + radius * radius);
    }
};

template<typename S, size_t Dim = 3>
class [[maybe_unused]] IsochroneSphere : public SphericalBase<IsochroneSphere<S, Dim>, S, Dim> {
    using Base = SphericalBase<IsochroneSphere<S, Dim>, S, Dim>;

public:
    using Vector = typename Base::Vector;
    using Scalar = typename Base::Scalar;

    explicit IsochroneSphere(ODIN_CONST Scalar mu,
            ODIN_CONST Scalar                  radius,
            ODIN_CONST Vector                 &position = Vector::Zero())
        : Base(mu, radius, position) {}

    static Scalar static_potential_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        ODIN_CONST Scalar r = rel_position.norm();
        return -mu / (radius + std::sqrt(r * r + radius * radius));
    }

    static Vector static_acceleration_impl(
            ODIN_CONST Vector &rel_position, ODIN_CONST Scalar mu, ODIN_CONST Scalar radius) {
        ODIN_CONST Scalar r = rel_position.norm();
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

    for (int i = 1; i <= n; ++i) { result *= i; }

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
            Scalar factor = sqrt((2 - (m == 0 ? 1 : 0)) * (2 * l + 1) * factorial<Scalar>(l - m)
                                 / factorial<Scalar>(l + m));
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
            Scalar factor = sqrt(1 / ((2 - (m == 0 ? 1 : 0)) * (2 * l + 1))
                                 * factorial<Scalar>(l + m) / factorial<Scalar>(l - m));
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
                Scalar factor = sqrt(1.0 / ((2 - (m == 0 ? 1 : 0)) * (2 * l + 1))
                                     * factorial<Scalar>(l + m) / factorial<Scalar>(l - m));
                inputClm[l][m] *= factor;
                inputSlm[l][m] *= factor;
            }
        }
    }

    static void denormalize_in_place(Clm &inputClm, Slm &inputSlm) {
        for (int l = 0; l <= Degree; ++l) {
            for (int m = 0; m <= std::min(l, Order); ++m) {
                Scalar factor = sqrt((2 - (m == 0 ? 1 : 0)) * (2 * l + 1) * factorial<Scalar>(l - m)
                                     / factorial<Scalar>(l + m));
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
                Scalar factor = sqrt(1.0 / ((2 - (m == 0 ? 1 : 0)) * (2 * l + 1))
                                     * factorial<Scalar>(l + m) / factorial<Scalar>(l - m));
                inputClm[l][m] *= factor;
                inputSlm[l][m] *= factor;
            }
        }
    }

    static void denormalize_in_place(Clm &inputClm, Slm &inputSlm) {
        for (int l = 0; l < inputClm.size(); ++l) {
            for (int m = 0; m < inputClm[l].size(); ++m) {
                Scalar factor = sqrt((2 - (m == 0 ? 1 : 0)) * (2 * l + 1) * factorial<Scalar>(l - m)
                                     / factorial<Scalar>(l + m));
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

    static std::pair<Clm, Slm> add_dimensions(
            Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
        Clm C_lm = inputC_lm;
        Slm S_lm = inputS_lm;
        add_dimensions_in_place(mu, R, C_lm, S_lm);
        return {C_lm, S_lm};
    }

    std::pair<Clm, Slm> remove_dimensions(
            Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
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

    std::pair<Clm, Slm> add_dimensions(
            Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
        Clm C_lm = inputC_lm;
        Slm S_lm = inputS_lm;
        add_dimensions_in_place(mu, R, C_lm, S_lm);
        return {C_lm, S_lm};
    }

    std::pair<Clm, Slm> remove_dimensions(
            Scalar mu, Scalar R, const Clm &inputC_lm, const Slm &inputS_lm) {
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
class [[maybe_unused]] SphericalHarmonics
    : public SphericalBase<SphericalHarmonics<S, Degree, Order, Dim>, S, Dim> {
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

    explicit SphericalHarmonics(const Scalar mu,
            const Scalar                     radius,
            const Clm                        Cnm = CoefficientsHolder<Scalar, Degree, Order>::Clm(),
            const Slm                        Snm = CoefficientsHolder<Scalar, Degree, Order>::Slm(),
            const Vector                    &position = Vector::Zero())
        : Base(mu, radius, position),
          Clm_(Cnm),
          Slm_(Snm) {
        if constexpr (Degree == -1 && Order == -1) {
            int degree = Clm_.size() - 1;
            int order  = Clm_.front().size() - 1;
            Plm_ = AssociatedLegendrePolynomialsHolder<Scalar, Degree, Order>::init(degree, order);
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
                Scalar P_lm         = mlib::assoc_legendre(l, m, sin_phi);
                Scalar term = P_lm * (Clm_[l][m] * cos_m_lambda + Slm_[l][m] * sin_m_lambda);
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


    Vector static_acceleration_impl(
            const Vector &rel_position, const Scalar mu, const Scalar radius) {
        Scalar r = rel_position.norm();
        if (r < radius) {
            return HomogeneousSphere<Scalar, Dim>::static_acceleration_impl(
                    rel_position, mu, radius);
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
        //            Scalar P_ll = mlib::assoc_legendre(l, l, sin_phi);
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
            for (size_t m = 0; m <= l; ++m) { Plm_[l][m] = mlib::assoc_legendre(l, m, sin_phi); }
        }

        for (size_t l = 0; l < Clm_.size(); ++l) {
            for (size_t m = 0; m <= l; ++m) {
                Scalar cos_m_lambda = std::cos(m * lambda);
                Scalar sin_m_lambda = std::sin(m * lambda);
                Scalar common_term  = (Clm_[l][m] * cos_m_lambda + Slm_[l][m] * sin_m_lambda);
                acceleration[0] -= ((l + 1) * Plm_[l][m] * common_term)
                                 * a_r_ratio_pow_l;// radial component

                if (l == 0 && m == 0) { continue; }// only radial component for l=0, m=0
                Scalar x       = sin_phi;          // We make the change of variables x = sin(phi)
                Scalar dx_dphi = cos_phi;          // Derivative of sin(phi) w.r.t. phi

                // Compute the derivative of the associated Legendre function w.r.t. x
                Scalar dPlm_dx
                        = (1 - x * x) * ((l - m + 1) * x * Plm_[l][m] - (l + 1) * Plm_[l + 1][m]);

                // Then, by the chain rule, compute the derivative of P_{l, m}(sin(phi)) w.r.t. phi
                Scalar dPlm_dphi = dPlm_dx * dx_dphi;
                acceleration[1] += (dPlm_dphi * common_term) * a_r_ratio_pow_l;

                if (l == 1 && m == 0) { continue; }
                acceleration[2] += m
                                 * (Plm_[l][m] / cos_phi
                                         * (-Clm_[l][m] * sin_m_lambda + Slm_[l][m] * cos_m_lambda))
                                 * a_r_ratio_pow_l;
            }
            a_r_ratio_pow_l *= radius / r;
        }
        acceleration *= r2_inv;
        return spherical_to_cartesian(acceleration, lambda, phi);
    }
};

// Verify that it adheres to the GravitationalModel concept
static_assert(is_gravitational_model_v<PointMass<double>>,
        "HollowSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<HollowSphere<double>>,
        "HollowSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<HomogeneousSphere<double>>,
        "HomogeneousSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<PlummerSphere<double>>,
        "PlummerSphere does not comply with GravitationalModel.");
static_assert(is_gravitational_model_v<IsochroneSphere<double>>,
        "IsochroneSphere does not comply with GravitationalModel.");

}// namespace odin::gravity


#endif// SPHERICAL_HPP
