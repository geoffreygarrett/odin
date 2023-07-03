//
// Created by geoffrey on 6/6/23.
//

#include <Eigen/Dense>
#include <tuple>
#include <functional>

// CRTP base class
template<typename Derived, typename Scalar, int Dim>
struct Augmentation {
    using State = Eigen::Matrix<Scalar, Dim, 1>;
    using Covariance = Eigen::Matrix<Scalar, Dim, Dim>;

    std::pair<State, Covariance> data;

    // Call the state_transition function on the derived class
    State state_transition(const State& x, Scalar dt) {
        return static_cast<Derived*>(this)->state_transition_impl(x, dt);
    }

    // Call the process_covariance function on the derived class
    Covariance process_covariance(const State& x, Scalar dt) {
        return static_cast<Derived*>(this)->process_covariance_impl(x, dt);
    }

    // Call the measurement_model function on the derived class
    State measurement_model(const State& x) {
        return static_cast<Derived*>(this)->measurement_model_impl(x);
    }
};

template<typename Scalar, int x_dim, int v_dim, int n_dim>
struct NoiseAugmentation : Augmentation<NoiseAugmentation<Scalar, x_dim, v_dim, n_dim>, Scalar, x_dim> {
    using Base = Augmentation<NoiseAugmentation<Scalar, x_dim, v_dim, n_dim>, Scalar, x_dim>;
    using State = typename Base::State;
    using Covariance = typename Base::Covariance;

    State state_transition_impl(const State& x, Scalar dt) {
        // Implement the state transition logic here for NoiseAugmentation
    }

    Covariance process_covariance_impl(const State& x, Scalar dt) {
        // Implement the process covariance logic here for NoiseAugmentation
    }

    State measurement_model_impl(const State& x) {
        // Implement the measurement model logic here for NoiseAugmentation
    }
};

template<typename Scalar, int param_dim>
struct ParameterAugmentation : Augmentation<ParameterAugmentation<Scalar, param_dim>, Scalar, param_dim> {
    using Base = Augmentation<ParameterAugmentation<Scalar, param_dim>, Scalar, param_dim>;
    using State = typename Base::State;
    using Covariance = typename Base::Covariance;

    State state_transition_impl(const State& x, Scalar dt) {
        // Implement the state transition logic here for ParameterAugmentation
    }

    Covariance process_covariance_impl(const State& x, Scalar dt) {
        // Implement the process covariance logic here for ParameterAugmentation
    }

    State measurement_model_impl(const State& x) {
        // Implement the measurement model logic here for ParameterAugmentation
    }
};

template<typename Scalar, typename... Augmentations>
class AugmentationHandler;

// Base case for AugmentationHandler
template<typename Scalar>
class AugmentationHandler<Scalar> {
    // base behavior here
};

template<typename Scalar, typename... Augmentations>
class AugmentationHandler {
public:
    AugmentationHandler() = default;

    std::tuple<typename Augmentations::State...>
    state_transition(const std::tuple<typename Augmentations::State...> &x, double dt) {
        return std::make_tuple((Augmentations().state_transition(std::get<typename Augmentations::State>(x), dt))...);
    }

    std::tuple<typename Augmentations::Covariance...>
    process_covariance(const std::tuple<typename Augmentations::State...> &x, double dt) {
        return std::make_tuple((Augmentations().process_covariance(std::get<typename Augmentations::State>(x), dt))...);
    }

    std::tuple<typename Augmentations::State...>
    measurement_model(const std::tuple<typename Augmentations::State...> &x) {
        return std::make_tuple((Augmentations().measurement_model(std::get<typename Augmentations::State>(x)))...);
    }

    std::tuple<typename Augmentations::State...>
    construct_augmented_state() {
        return std::make_tuple(Augmentations().data.first...);
    }

    std::tuple<typename Augmentations::Covariance...>
    construct_augmented_covariance() {
        return std::make_tuple(Augmentations().data.second...);
    }
};