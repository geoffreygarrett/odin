#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include "Eigen/Dense"
#include <array>
#include <functional>
#include <odin/core/numerics/integrators/base_integrator.hpp>
#include <odin/core/numerics/integrators/butcher_tableau.hpp>
#include <odin/core/policy/cache.hpp>
#include <tuple>
#include <vector>


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

    static constexpr auto a           = Tableau::a;          // Runge-Kutta matrix
    static constexpr auto b           = Tableau::b;          // weights
    static constexpr auto c           = Tableau::c;          // nodes
    static constexpr auto stages      = Tableau::stages;     // number of stages
    static constexpr auto is_embedded = Tableau::is_extended;// whether the tableau is extended

    CachePolicy cache_policy;

    /**
     * @brief Constructor that takes the dynamics function as a parameter.
     *
     * @param dynamics_function The dynamics function.
     */
    // Template constructor for efficiency
    //    template<typename Callable>
    //    explicit RungeKuttaIntegrator(Callable &&dynamics_function)
    //        : Base(std::forward<Callable>(dynamics_function)) {}

    explicit RungeKuttaIntegrator(const Function &dynamics_function)
        : Base(dynamics_function) {}

    // This function will only be considered if Tableau::is_extended, therefore rk can be embedded
    template<typename T = Tableau>
    State step_impl(const State &state, Scalar dt) {
        if (!cache_policy.compare_and_store_if_invalid_impl(state, dt)) {
            calculate_k(state, dt);
        }
        State new_state = state;
        step_helper(new_state, dt, std::make_index_sequence<stages>{});

        if constexpr (Tableau::is_extended) {
            State z_next = state;
            adjust_step_size_impl_helper(z_next, dt, state, std::make_index_sequence<stages>{});
            z_next_ = z_next;
        }

        return new_state;
    }


    template<typename T = Tableau>
    std::pair<State, Scalar> step_with_dt(const State &state, Scalar dt) {
        // Implement the actual step using the step_impl function
        State new_state = step_impl(state, dt);

        // Adjust step size, if the tableau is adaptive
        Scalar new_dt = dt;
        if constexpr (Tableau::is_extended) {
            new_dt = adjust_step_size_impl(dt, state, new_state);
        }

        // Return the new state and the new step size
        return {new_state, new_dt};
    }

    // https://www.researchgate.net/publication/346650178_Runge_Kutta_Fehlberg_45_MATLAB_code
    State z_next_;// To store z_next calculated in step_impl


    template<std::size_t... Is>
    void step_helper(State &new_state, Scalar dt, std::index_sequence<Is...>) {
        ((new_state += k_cache[Is] * c[Is] * dt), ...);
    }

    template<std::size_t... Is>
    void step_impl_helper(State &new_state, Scalar dt, std::index_sequence<Is...>) {
        static constexpr auto a = Tableau::a;
        static constexpr auto b = Tableau::b;
        if constexpr (Tableau::is_extended) {
            static constexpr auto b_star = Tableau::b_star;
            ((new_state += dt * (b[Is] * k_cache[Is] + b_star[Is] * k_cache[Is])), ...);
        } else {
            ((new_state += dt * b[Is] * k_cache[Is]), ...);
        }
    }


    template<std::size_t... Is>
    void adjust_step_size_impl_helper(State &z_next, Scalar dt, const State &state,
                                      std::index_sequence<Is...>) {
        static constexpr auto b_star = Tableau::b_star;
        ((z_next += dt * b_star[Is] * k_cache[Is]), ...);
    }

    // This function will only be considered if Tableau::is_embedded is true
    template<typename T = Tableau>
    std::enable_if_t<T::is_embedded, Scalar>
    adjust_step_size_impl(Scalar dt, const State &state, Scalar /*current_time*/) {
        const Scalar tol       = 1e-6;// Set the desired tolerance
        const Scalar abs_error = (z_next_ - state).norm();
        const Scalar s         = 0.84 * std::pow(tol * dt / abs_error, 0.25);
        // TODO: Check if 0.84 should be used or the analytical value it approximates.
        return dt * s;// The optimal step size
    }

    template<std::size_t... Is>
    void calculate_error_helper(Error &error, std::index_sequence<Is...>) {
        static constexpr auto b_star = Tableau::b_star;
        ((error += (c[Is] - b_star[Is]) * dt_cache * k_cache[Is]), ...);
    }

    // This function will only be considered if Tableau::is_embedded is true
    template<typename T = Tableau>
    std::enable_if_t<T::is_embedded, Error>
    estimate_error_impl(const State &initial, Scalar dt) {
        if (!cache_policy.compare_and_store_if_invalid_impl(initial, dt)) {
            calculate_k(initial, dt);
        }
        Error error = Error::Zero(initial.size());
        calculate_error_helper<stages>(error, std::make_index_sequence<stages>{});
        return error;
    }

    // Setter for tolerance
    void set_tolerance(Scalar new_tol) {
        tolerance = new_tol;
    }

    // Getter for tolerance
    Scalar get_tolerance() const {
        return tolerance;
    }

private:
    Scalar tolerance;// tolerance for error control

    template<std::size_t j, std::size_t i, std::size_t... Js>
    State calculate_state_derivative(
            const State &initial,
            Scalar       dt,
            State        state_derivative,
            std::index_sequence<Js...>) {
        ((state_derivative += k_cache[Js] * a[i][Js]), ...);
        return state_derivative;
    }

    template<std::size_t i>
    void calculate_k_iteration(const State &initial, Scalar dt) {
        if constexpr (i != 0) {
            State state_derivative = State::Zero(initial.size());
            state_derivative       = calculate_state_derivative<0, i>(
                    initial,
                    dt,
                    state_derivative,
                    std::make_index_sequence<i>{});
            k_cache[i] = this->dynamics_function_(initial + dt * state_derivative, c[i] * dt);
        }
    }

    template<std::size_t... Is>
    void calculate_k_inner(const State &initial, Scalar dt, std::index_sequence<Is...>) {
        ((calculate_k_iteration<Is>(initial, dt)), ...);
    }

    void calculate_k(const State &initial, Scalar dt) {
        cache_is_valid_ = true;
        k_cache[0]      = this->dynamics_function_(initial, Scalar(0));
        calculate_k_inner(initial, dt, std::make_index_sequence<stages>{});
        state_cache = initial;
        dt_cache    = dt;
    }

    ButcherTableau<Scalar, Tableau::stages> tableau;
    std::array<Derivative, stages>          k_cache;
    State                                   state_cache;
    Scalar                                  dt_cache;
    bool                                    cache_is_valid_ = false;
};

//template<typename State,
//         typename Function,
//         typename Tableau,
//         typename CachePolicy = NoCachePolicy<State, typename State::Scalar>>
//class RungeKuttaIntegrator
//    : public BaseIntegrator<
//              RungeKuttaIntegrator<State, Function, Tableau, CachePolicy>, State, Function> {
//public:
//    using Base = BaseIntegrator<RungeKuttaIntegrator<State, Function, Tableau, CachePolicy>, State, Function>;
//    using Base::Base;
//
//    using typename Base::Derivative;
//    using typename Base::Error;
//    using typename Base::Scalar;
//
//    static constexpr auto rk_matrix   = Tableau::a;          // Runge-Kutta matrix
//    static constexpr auto weights     = Tableau::b;          // weights
//    static constexpr auto nodes       = Tableau::c;          // nodes
//    static constexpr auto stages      = Tableau::stages;     // number of stages
//    static constexpr auto is_embedded = Tableau::is_extended;// whether the tableau is extended
//
//    CachePolicy cache_policy;
//
//    /**
//     * @brief Constructor that takes the dynamics function as a parameter.
//     *
//     * @param dynamics_function The dynamics function.
//     */
//    // Template constructor for efficiency
//    //    template<typename Callable>
//    //    explicit RungeKuttaIntegrator(Callable &&dynamics_function)
//    //        : Base(std::forward<Callable>(dynamics_function)) {}
//
//    explicit RungeKuttaIntegrator(const Function &derivative_function)
//        : Base(derivative_function) {}
//
//    // This function will only be considered if Tableau::is_extended, therefore rk can be embedded
//    template<typename T = Tableau>
//    State step_impl(const State &state, Scalar dt) {
//        if (!cache_policy.compare_and_store_if_invalid_impl(state, dt)) {
//            calculate_k(state, dt);
//        }
//        State new_state = state;
//        step_helper(new_state, dt, std::make_index_sequence<stages>{});
//
//        if constexpr (Tableau::is_extended) {
//            State z_next = state;
//            adjust_step_size_impl_helper(z_next, dt, state, std::make_index_sequence<stages>{});
//            z_next_ = z_next;
//        }
//
//        return new_state;
//    }
//
//
//    template<typename T = Tableau>
//    std::pair<State, Scalar> step_with_dt(const State &state, Scalar dt) {
//        // Implement the actual step using the step_impl function
//        State new_state = step_impl(state, dt);
//
//        // Adjust step size, if the tableau is adaptive
//        Scalar new_dt = dt;
//        if constexpr (Tableau::is_extended) {
//            new_dt = adjust_step_size_impl(dt, state, new_state);
//        }
//
//        // Return the new state and the new step size
//        return {new_state, new_dt};
//    }
//
//    // https://www.researchgate.net/publication/346650178_Runge_Kutta_Fehlberg_45_MATLAB_code
//    State z_next_;// To store z_next calculated in step_impl
//
//
//    template<std::size_t... Is>
//    void step_helper(State &new_state, Scalar dt, std::index_sequence<Is...>) {
//        ((new_state += k_cache[Is] * c[Is] * dt), ...);
//    }
//
//    template<std::size_t... Is>
//    void step_impl_helper(State &new_state, Scalar dt, std::index_sequence<Is...>) {
//        static constexpr auto a = Tableau::a;
//        static constexpr auto b = Tableau::b;
//        if constexpr (Tableau::is_extended) {
//            static constexpr auto b_star = Tableau::b_star;
//            ((new_state += dt * (b[Is] * k_cache[Is] + b_star[Is] * k_cache[Is])), ...);
//        } else {
//            ((new_state += dt * b[Is] * k_cache[Is]), ...);
//        }
//    }
//
//
//    template<std::size_t... Is>
//    void adjust_step_size_impl_helper(State &z_next, Scalar dt, const State &state,
//                                      std::index_sequence<Is...>) {
//        static constexpr auto b_star = Tableau::b_star;
//        ((z_next += dt * b_star[Is] * k_cache[Is]), ...);
//    }
//
//    // This function will only be considered if Tableau::is_embedded is true
//    template<typename T = Tableau>
//    std::enable_if_t<T::is_embedded, Scalar>
//    adjust_step_size_impl(Scalar dt, const State &state, Scalar /*current_time*/) {
//        const Scalar tol       = 1e-6;// Set the desired tolerance
//        const Scalar abs_error = (z_next_ - state).norm();
//        const Scalar s         = 0.84 * std::pow(tol * dt / abs_error, 0.25);
//        // TODO: Check if 0.84 should be used or the analytical value it approximates.
//        return dt * s;// The optimal step size
//    }
//
//    template<std::size_t... Is>
//    void calculate_error_helper(Error &error, std::index_sequence<Is...>) {
//        static constexpr auto b_star = Tableau::b_star;
//        ((error += (c[Is] - b_star[Is]) * dt_cache * k_cache[Is]), ...);
//    }
//
//    // This function will only be considered if Tableau::is_embedded is true
//    template<typename T = Tableau>
//    std::enable_if_t<T::is_embedded, Error>
//    estimate_error_impl(const State &initial, Scalar dt) {
//        if (!cache_policy.compare_and_store_if_invalid_impl(initial, dt)) {
//            calculate_k(initial, dt);
//        }
//        Error error = Error::Zero(initial.size());
//        calculate_error_helper<stages>(error, std::make_index_sequence<stages>{});
//        return error;
//    }
//
//    // Setter for tolerance
//    void set_tolerance(Scalar new_tol) {
//        tolerance = new_tol;
//    }
//
//    // Getter for tolerance
//    Scalar get_tolerance() const {
//        return tolerance;
//    }
//
//private:
//    Scalar tolerance;// tolerance for error control
//
//    template<std::size_t j, std::size_t i, std::size_t... Js>
//    State calculate_state_derivative(
//            const State &initial,
//            Scalar       dt,
//            State        state_derivative,
//            std::index_sequence<Js...>) {
//        ((state_derivative += k_cache[Js] * rk_matrix[i][Js]), ...);
//        return state_derivative;
//    }
//
//    template<std::size_t i>
//    void calculate_k_iteration(const State &initial, Scalar dt) {
//        if constexpr (i != 0) {
//            State state_derivative = State::Zero(initial.size());
//            state_derivative       = calculate_state_derivative<0, i>(
//                    initial,
//                    dt,
//                    state_derivative,
//                    std::make_index_sequence<i>{});
//            k_cache[i] = this->dynamics_function_(initial + dt * state_derivative, nodes[i] * dt);
//        }
//    }
//
//    template<std::size_t... Is>
//    void calculate_k_inner(const State &initial, Scalar dt, std::index_sequence<Is...>) {
//        ((calculate_k_iteration<Is>(initial, dt)), ...);
//    }
//
//    void calculate_k(const State &initial, Scalar dt) {
//        cache_is_valid_ = true;
//        k_cache[0]      = this->dynamics_function_(initial, Scalar(0));
//        calculate_k_inner(initial, dt, std::make_index_sequence<stages>{});
//        state_cache = initial;
//        dt_cache    = dt;
//    }
//
//    ButcherTableau<Scalar, Tableau::stages> tableau;
//    std::array<Derivative, stages>          k_cache;
//    State                                   state_cache;
//    Scalar                                  dt_cache;
//    bool                                    cache_is_valid_ = false;
//};


#endif// RUNGE_KUTTA_HPP
