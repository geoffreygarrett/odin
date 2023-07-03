#ifndef OBSERVATION_MODEL_HPP
#define OBSERVATION_MODEL_HPP

#include "Eigen/Core"
#include "autodiff/forward/real/eigen.hpp"

template<typename Derived, typename StateType, typename MeasurementType, typename JacobianType, typename ParamType>
struct MeasurementModel {

    static constexpr int state_dim = StateType::RowsAtCompileTime;
    static constexpr int param_dim = ParamType::RowsAtCompileTime;
    static constexpr int meas_dim = MeasurementType::RowsAtCompileTime;

    constexpr MeasurementType compute_measurement(const StateType &state, const ParamType &param) const {
        return static_cast<const Derived *>(this)->compute_measurement_impl(state, param);
    }

    template<typename T = Derived>
    constexpr auto compute_jacobian(const StateType &state, const ParamType &params) const
    -> std::enable_if_t<has_compute_analytical_jacobian<T, StateType, JacobianType>(), JacobianType> {
        return static_cast<const T *>(this)->compute_analytical_jacobian(state, params);
    }

    template<typename T = Derived>
    constexpr auto compute_jacobian(const StateType &state, const ParamType &params) const
    -> std::enable_if_t<!has_compute_analytical_jacobian<T, StateType, JacobianType>(), JacobianType> {
        return compute_autodiff_jacobian(state, params);
    }
    

    constexpr JacobianType compute_autodiff_jacobian(const StateType &state, const ParamType &params) const {
        using namespace autodiff;
        using Scalar = typename MeasurementType::Scalar;

        using ADStateType = Eigen::Matrix<autodiff::real, state_dim, 1>;
        ADStateType state_ad = state.template cast<autodiff::real>();

        using ADParamType = Eigen::Matrix<autodiff::real, param_dim, 1>;
        ADParamType params_ad = params.template cast<autodiff::real>();

        // Wrapper lambda function to pass into autodiff::jacobian
        auto f = [&](const ADStateType &state, const ADParamType &params) {
            return static_cast<const Derived *>(this)->compute_measurement_impl(state, params);
        };

        MeasurementType result;   // the output vector F = f(state, params) evaluated together with Jacobian below

        // Compute the Jacobian
        JacobianType jacobian = autodiff::jacobian(f, wrt(state_ad, params_ad), at(state_ad, params_ad), result);

        return jacobian;
    }
};

#endif // OBSERVATION_MODEL_HPP
