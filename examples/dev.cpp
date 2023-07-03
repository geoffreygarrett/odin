#include "include/odin/core/numerics/integrators/base_integrator.hpp"
#include "include/odin/core/policy/cache.hpp"
#include <array>
#include <deque>
#include <iostream>
#include <odin/logging.hpp>

template<typename Scalar, size_t Stages, typename ExtensionType = void>
struct ButcherTableau {
    static constexpr size_t stages      = Stages;
    static constexpr bool   is_extended = false;

    std::array<std::array<Scalar, Stages>, Stages> a;// Runge-Kutta matrix
    std::array<Scalar, Stages>                     b;// weights
    std::array<Scalar, Stages>                     c;// nodes
};

template<typename Scalar, size_t Stages>
struct ButcherTableau<Scalar, Stages, std::array<Scalar, Stages>> {
    static constexpr size_t stages      = Stages;
    static constexpr bool   is_extended = true;

    std::array<std::array<Scalar, Stages - 1>, Stages> a;     // Runge-Kutta matrix
    std::array<Scalar, Stages>                         b;     // weights
    std::array<Scalar, Stages>                         c;     // nodes
    std::array<Scalar, Stages>                         b_star;// lower-order weights
};

namespace butcher_tableau {
    // A non-adaptive tableau
    // Forward Euler method
    template<typename T = double>
    struct Euler : public ButcherTableau<T, 1> {
        constexpr static std::array<std::array<T, 1>, 1> a = {{{0.0}}};
        constexpr static std::array<T, 1>                c = {0.0};
        constexpr static std::array<T, 1>                b = {1.0};
    };

    template<typename T = double>
    struct Heun : public ButcherTableau<T, 2> {
        constexpr static std::array<std::array<T, 2>, 2> a = {
                {{0.0, 0.0}, {1.0, 0.0}}
        };
        constexpr static std::array<T, 2> c = {0.0, 1.0};
        constexpr static std::array<T, 2> b = {0.5, 0.5};
    };

    // Heun's method with Euler method as its lower order method
    template<typename T = double>
    struct HeunEuler : public ButcherTableau<T, 2, std::array<T, 2>> {
        constexpr static auto             a      = Heun<T>::a;
        constexpr static auto             c      = Heun<T>::c;
        constexpr static auto             b      = Heun<T>::b;
        constexpr static std::array<T, 2> b_star = {1.0, 0.0};
    };


    template<typename T = double>
    struct Fehlberg : public ButcherTableau<T, 6> {
        constexpr static std::array<std::array<T, 6>, 6> a = {
                {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {1.0 / 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0, 0.0},
                 {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0, 0.0},
                 {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0, 0.0},
                 {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0, 0.0}}
        };
        constexpr static std::array<T, 6> c = {
                0.0,
                1.0 / 4.0,
                3.0 / 8.0,
                12.0 / 13.0,
                1.0,
                1.0 / 2.0};
        constexpr static std::array<T, 6> b = {
                16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};
        constexpr static std::array<T, 6> b_star = {
                25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0};
    };


}// namespace butcher_tableau

template<typename T>
concept Scalar = requires(T v) {
    { std::abs(v) } -> std::same_as<T>;
};

template<typename T>
concept Vector = requires(T v) {
    { v.norm() } -> std::same_as<T>;
};

template<typename T>
concept Tolerance = requires(T t, typename T::value_type error, typename T::value_type value) {
    { t.within_tolerance(error, value) } -> std::same_as<bool>;
};

template<Scalar S>
class AbsoluteTolerance {
public:
    using value_type = S;

    bool within_tolerance(const S &error, const S &value) const {
        return std::abs(error) <= value;
    }

    template<Vector V>
    bool within_tolerance(const V &error, const S &value) const {
        return error.norm() <= value;
    }
};

template<Scalar S>
class RelativeTolerance {
public:
    using value_type = S;

    bool within_tolerance(const S &error, const S &value) const {
        return std::abs(error / value) <= 1;
    }

    template<Vector V>
    bool within_tolerance(const V &error, const S &value) const {
        return (error.norm() / value) <= 1;
    }
};

template<Scalar S>
class DefaultTolerance {
public:
    using value_type = S;

    bool within_tolerance(const S &error, const S &value) const {
        return std::abs(error) <= (std::numeric_limits<S>::epsilon() * std::abs(value));
    }

    template<Vector V>
    bool within_tolerance(const V &error, const S &value) const {
        return error.norm() <= (std::numeric_limits<S>::epsilon() * std::abs(value));
    }
};

/**
 * @brief Class that implements a Runge-Kutta ODE solver with a configurable Butcher tableau and event handling mechanism
 *
 * @tparam State The type of the state vector
 * @tparam Derivative The type of the derivative vector
 * @tparam N The number of stages in the Runge-Kutta method
 * @tparam EventPolicy The type of the event handling mechanism
 * @tparam Tableau The Butcher tableau of the Runge-Kutta method
 */
template<typename State,
         typename Function,
         typename Tableau,
         typename CachePolicy = NoCachePolicy<State, typename State::Scalar>>
class RungeKuttaIntegrator
    : public BaseIntegrator<
              RungeKuttaIntegrator<State, Function, Tableau, CachePolicy>, State, Function> {
public:
    using Base = BaseIntegrator<RungeKuttaIntegrator<State, Function, Tableau, CachePolicy>, State, Function>;
    using Base::Base;

    using typename Base::Derivative;
    using typename Base::Error;
    using typename Base::Scalar;

    static constexpr auto rk_matrix   = Tableau::a;          // Runge-Kutta matrix
    static constexpr auto rk_weights  = Tableau::b;          // weights
    static constexpr auto rk_nodes    = Tableau::c;          // nodes
    static constexpr auto rk_stages   = Tableau::stages;     // number of stages
    static constexpr auto is_embedded = Tableau::is_extended;// whether the tableau is extended

    CachePolicy cache_policy;

    // Common to all ordinary differential equation methods ////////////////////////////////////////
    template<typename T = Tableau>
    std::enable_if_t<!T::is_extended, State>
    step_impl(
            const State &state,
            Scalar       t0,
            Scalar       tf,
            Scalar       dt) {
        return primary_step(state, t0, tf, dt);
    }

    template<typename T = Tableau>
    std::enable_if_t<T::is_extended, State>
    step_impl(
            const State &state,
            Scalar       t0,
            Scalar       tf,
            Scalar       dt) {
        return adaptive_step(state, t0, tf, dt);
    }

    State integrate_impl(const State &initial_state, Scalar t0, Scalar tf, Scalar dt) {
        // Initialize the state and time
        State state = initial_state;
        t_          = t0;

        // Initialize the cache
        cache_policy.initialize(initial_state);

        // Integrate until the final time is reached
        while (t_ < tf) {
            // Compute the derivatives at all stages
            compute_all_stage_derivatives(state, dt);

            // Compute the next state
            state = compute_next_state(state, dt);

            // Update the time
            t_ += dt;
        }

        // Return the final state
        return state;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    template<std::size_t... rows>
    void step_helper(State &new_state, Scalar dt, std::index_sequence<rows...>) {
        ((new_state += stage_derivatives_[rows] * rk_nodes[rows] * dt), ...);
    }

    template<std::size_t... rows>
    State aggregate_weights(const std::array<Scalar, rk_stages> &weights, Scalar dt, std::index_sequence<rows...>) {
        return (accumulate_weight<rows>(weights, dt) + ...);
    }

    template<std::size_t row>
    State accumulate_weight(const std::array<Scalar, rk_stages> &weights, Scalar dt) {
        if constexpr (weights[row] != 0) {
            return stage_derivatives_[row] * weights[row] * dt;
        } else {
            return State::Zero();
            // TODO: Could be optimised (memory-wise, maybe through a sentinel value)
        }
    }

    std::array<Derivative, rk_stages> get_stage_derivatives() const {
        return stage_derivatives_;
    }

    State
    primary_step(
            const State &state,
            Scalar       t0,
            Scalar       tf,
            Scalar       dt) {
        compute_all_stage_derivatives(state, dt);
        State new_state = state;
        return aggregate_weights<rk_stages>(rk_weights, dt, std::make_index_sequence<rk_stages>{});
    }

    template<typename T = Tableau>
    std::enable_if_t<T::is_extended, State>
    refined_step(
            const State &state,
            Scalar       t0,
            Scalar       tf,
            Scalar       dt) {
        compute_all_stage_derivatives(state, dt);
        State new_state = state;
        return aggregate_weights<rk_stages>(T::b_star, dt, std::make_index_sequence<rk_stages>{});
    }


    template<typename T = Tableau>
    std::enable_if_t<T::is_extended, Scalar> const compute_scale(const State &error) {
        return pow((tolerance_ / (2 * error.norm())), 1 / 4);
    }

    template<typename T = Tableau>
    std::enable_if_t<T::is_extended, std::tuple<State, Scalar>>
    adaptive_step(
            const State &state,
            Scalar       t0,
            Scalar       tf,
            Scalar       dt,
            Scalar       tolerance     = 1e-6,
            unsigned int max_iter      = 1000,
            Scalar       min_step_size = 1e-6
    ) {
        unsigned int iteration = 0;

        while (true) {
            compute_all_stage_derivatives(state, dt);
            auto primary = primary_step(state, t0, tf, dt);
            auto error   = refined_step(state, t0, tf, dt) - primary;
            auto new_dt  = compute_scale(error) * dt;

            if (++iteration >= max_iter) {
                ODIN_LOG_WARNING << "Maximum number of iterations exceeded in adaptive Runge-Kutta, returning current state and time step";
                return std::make_tuple(primary, dt);
            }

            if (new_dt < min_step_size) {
                ODIN_LOG_WARNING << "Minimum step size reached in adaptive Runge-Kutta, returning current state and time step";
                return std::make_tuple(primary, dt);
            }

            if (error.norm() < tolerance_) {
                return std::make_tuple(primary, new_dt);
            }

            // Update the time step for the next iteration.
            dt = new_dt;
        }
    }

    /**
     * @brief Integrates the ODE from the initial state to the final state
     *
     * @param initial_state The initial state
     * @param t0 The initial time
     * @param tf The final time
     * @param dt The time step
     * @return std::pair<State, Error> The final state and the error
     */
    State integrate(const State &initial_state, Scalar t0, Scalar tf, Scalar dt) {
        // Initialize the state and time
        State state = initial_state;
        t_          = t0;

        // Initialize the cache
        cache_policy.initialize(initial_state);

        // Integrate until the final time is reached
        while (t_ < tf) {
            // Compute the derivatives at all stages
            compute_all_stage_derivatives(state, dt);

            // Compute the next state
            state = compute_next_state(state, dt);

            // Update the time
            t_ += dt;

            // Update the cache
            cache_policy.update(state);
        }

        // Return the final state and the error
        return state;
    }


private:
    Scalar tolerance_;
    // k1, k2, k3, ...
    void compute_all_stage_derivatives(const State &initial_state, Scalar dt) {
        stage_derivatives_[0] = this->differential_equation_(initial_state, t_);
        compute_stage_derivatives(initial_state, dt, std::make_index_sequence<rk_stages>{});
    }

    template<std::size_t... rows>
    void compute_stage_derivatives(// k1, k2, k3, ...
            const State &initial_state,
            Scalar       dt,
            std::index_sequence<rows...>) {
        ((calculate_stage_derivative<rows>(initial_state, dt)), ...);
    }

    template<std::size_t row>
    void calculate_stage_derivative(const State &initial, Scalar dt) {//ki
        if constexpr (row != 0) {
            State derivative = State::Zero(initial.size());
            derivative       = calculate_derivative<0, row>(
                    initial,
                    dt,
                    derivative,
                    std::make_index_sequence<row>{});
            stage_derivatives_[row] = this->dynamics_function_(
                    initial + dt * derivative, t_ + rk_nodes[row] * dt);
        }
    }

    template<std::size_t j, std::size_t row, std::size_t... cols>
    State calculate_derivative(
            const State &initial,
            Scalar       dt,
            State        derivative,
            std::index_sequence<cols...>) {
        ((derivative += stage_derivatives_[cols] * rk_matrix[row][cols]), ...);
        return derivative;
    }

    ButcherTableau<Scalar, rk_stages> tableau_;            // The Butcher tableau
    std::array<Derivative, rk_stages> stage_derivatives_{};// k1, k2, k3, ...
    State                             state_;              // The current state
    Scalar                            dt_;                 // The current step size
    Scalar                            t_;                  // The current independent variable
};
