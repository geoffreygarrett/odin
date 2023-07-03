#ifndef OBSERVATION_MODEL_HPP
#define OBSERVATION_MODEL_HPP

#include "Eigen/Core"
#include "autodiff/forward/dual.hpp"
#include <type_traits>

// C++20: Concepts check
// ---------------------
// C++17 replacement would be a std::enable_if in combination with std::is_same and type traits
// It would look like this:
// template<typename T, typename = void>
// struct has_compute_analytical_jacobian : std::false_type {};
//
// template<typename T>
// struct has_compute_analytical_jacobian<T, std::void_t<decltype(std::declval<T>().compute_analytical_jacobian(
//         std::declval<X>()))>> : std::true_type {};
//
// template<typename T>
// inline constexpr bool has_compute_analytical_jacobian_v = has_compute_analytical_jacobian<T>::value;
template<typename T, typename X, typename JacobianType>
concept has_compute_analytical_jacobian = requires(T a, X x) {
    { a.compute_analytical_jacobian(x) } -> std::same_as<JacobianType>;
};

// Requires C++17 or later
template<typename Derived, typename X, typename Z, typename Q, int JRows = Z::RowsAtCompileTime, int JCols = X::RowsAtCompileTime>
struct MeasurementModel {
    // Using C++17 type alias to improve code readability
    using JacobianType = Eigen::Matrix<double, JRows, JCols>;

    // C++14: Introduced variable templates
    static constexpr int StateDimension = X::RowsAtCompileTime;
    static constexpr int MeasurementDimension = Z::RowsAtCompileTime;

    // Measurement Model: z = h(x, q) + v, where:
    // z is the measurement
    // x is the state
    // q are the parameters affecting the measurement
    // h() is the function mapping the state and parameters to the measurement
    // v is the measurement noise
    // C++14: Constexpr member functions allow computations to be performed at compile-time, improving performance
    constexpr Z compute_measurement(const X &state, const Q &params) const {
        return static_cast<const Derived *>(this)->compute_measurement_impl(state, params);
    }

    // C++14: Constexpr member functions allow computations to be performed at compile-time, improving performance
    constexpr JacobianType compute_jacobian(const X &state, const Q &params, double delta = 1e-8) const {
        if constexpr (has_compute_analytical_jacobian<Derived, X, JacobianType>) {
            return static_cast<const Derived *>(this)->compute_analytical_jacobian(state, params);
        } else {
            return compute_numerical_jacobian(state, params, delta);
        }
    }
    // C++14: Constexpr member functions allow computations to be performed at compile-time, improving performance
    constexpr JacobianType compute_numerical_jacobian(const X &state, const Q &params, double delta = 1e-8) const {
        JacobianType jacobian;
        X perturbation = X::Zero();

        for (int i = 0; i < state.size(); ++i) {
            perturbation[i] = delta;
            const Z measurement_plus = compute_measurement(state + perturbation, params);
            const Z measurement_minus = compute_measurement(state - perturbation, params);

            jacobian.col(i) = (measurement_plus - measurement_minus) / (2 * delta);

            perturbation[i] = 0.0;
        }

        return jacobian;
    }

    // The Jacobian calculation using autodiff
    template<typename F>
    constexpr JacobianType compute_autodiff_jacobian(const F& f, const X &state, const Q &params) const {
        JacobianType jacobian;

        for (int i = 0; i < state.size(); ++i) {
            autodiff::dual x = state[i];
            autodiff::dual u = f(x, params); // assuming your measurement function can work with autodiff dual type

            jacobian.col(i) = autodiff::derivative(f, wrt(x), at(x, params));
        }

        return jacobian;
    }

    // Providing a way of verifying if the Jacobian is numerical/analytical via compile time trait
    template<typename T = Derived>
    constexpr bool is_analytical_jacobian() {
        return has_compute_analytical_jacobian<T, X, JacobianType>;
    }

    template<typename T = Derived>
    constexpr bool is_numerical_jacobian() {
        return !has_compute_analytical_jacobian<T, X, JacobianType>;
    }


};

#endif // OBSERVATION_MODEL_HPP