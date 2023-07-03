#ifndef OBSERVATION_MODEL_HPP
#define OBSERVATION_MODEL_HPP

#include "Eigen/Core"
#include "autodiff/forward/real/eigen.hpp"

template<typename T, typename X, typename JacobianType>
concept has_compute_analytical_jacobian = requires(T a, X x) {
    { a.compute_analytical_jacobian(x) } -> std::same_as<JacobianType>;
};

template<typename Derived, typename StateType, typename MeasurementType, typename JacobianType>
struct MeasurementModel {

//    static constexpr int state_dim = StateType::RowsAtCompileTime;
//    static constexpr int meas_dim = MeasurementType::RowsAtCompileTime;

    constexpr MeasurementType compute_measurement(const StateType &state) const {
        return static_cast<const Derived *>(this)->compute_measurement_impl(state);
    }

    template<typename T = Derived>
    constexpr auto compute_jacobian(const StateType &state) const
    -> std::enable_if_t<has_compute_analytical_jacobian<T, StateType, JacobianType>, JacobianType> {
        return static_cast<const T *>(this)->compute_analytical_jacobian(state);
    }

    template<typename T = Derived>
    constexpr auto compute_jacobian(const StateType &state) const
    -> std::enable_if_t<!has_compute_analytical_jacobian<T, StateType, JacobianType>, JacobianType> {
        return compute_autodiff_jacobian(state);
    }


    constexpr JacobianType compute_autodiff_jacobian(const StateType &state) const {
        using namespace autodiff;
        using Scalar = typename MeasurementType::Scalar;

        using ADStateType = Eigen::Matrix<autodiff::real, state_dim, 1>;
        ADStateType state_ad = state.template cast<autodiff::real>();

        // Wrapper lambda function to pass into autodiff::jacobian
        auto f = [&](const ADStateType &state) {
            return static_cast<const Derived *>(this)->compute_measurement_impl(state);
        };

        MeasurementType result;   // the output vector F = f(state, params) evaluated together with Jacobian below

        // Compute the Jacobian
        JacobianType jacobian = autodiff::jacobian(f, wrt(state_ad), at(state_ad), result);

        return jacobian;
    }
};

#endif // OBSERVATION_MODEL_HPP
