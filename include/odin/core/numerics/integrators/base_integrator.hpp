#ifndef BASE_INTEGRATOR_HPP
#define BASE_INTEGRATOR_HPP

#include <deque>
#include <odin/state/state.hpp>


// Base template
template<typename State, typename Enable = void>
struct ScalarType {
    using type = State;
};

// Specialization for types that have a Scalar typedef
template<typename State>
struct ScalarType<State, typename std::enable_if_t<std::is_class_v<State>>> {
    using type = typename State::Scalar;
};

template<typename Derived, typename State, typename Function>
class BaseIntegrator {
public:
    using Derivative = State;
    using Scalar     = typename ScalarType<State>::type;
    using Error      = State;

    BaseIntegrator() = default;

    template<typename Callable>
    explicit BaseIntegrator(Callable &&dynamics_function)
        : dynamics_function_(std::forward<Callable>(dynamics_function)) {
//        static_assert(std::is_same_v<State, std::invoke_result_t<Callable, State, Scalar>>,
//                      "The Callable must return the same type as the State");
    }

    template<typename Func>
    void set_dynamics_function(Func &&dynamics_function) {
        dynamics_function_ = std::forward<Func>(dynamics_function);
//        static_assert(std::is_same_v<State, std::invoke_result_t<Func, State, Scalar>>,
//                      "The Callable must return the same type as the State");
    }

    Scalar adjust_step_size(Scalar dt, const State &state, Scalar current_time) {
        if constexpr (requires { static_cast<Derived *>(this)->adjust_step_size_impl(dt, state, current_time); }) {
            // Call the derived class' implementation
            return static_cast<Derived *>(this)->adjust_step_size_impl(dt, state, current_time);
        } else {
            // By default, don't change the step size
            return dt;
        }
    }

    [[maybe_unused]] Scalar estimate_error(const State &state, Scalar dt) {
        if constexpr (requires { static_cast<Derived *>(this)->estimate_error_impl(state, dt); }) {
            // Call the derived class' implementation
            return static_cast<Derived *>(this)->estimate_error_impl(state, dt);
        } else {
            // By default, report no error
            return Scalar(0);
        }
    }

    State integrate(State &state, Scalar dt, Scalar t_final, Function dynamics_function) {
        set_dynamics_function(dynamics_function);
        return integrate(state, dt, t_final);
    }

    State integrate(State &state, Scalar dt, Scalar t_final, Scalar t_initial = Scalar(0)) {
        Scalar t_current = t_initial;
        while (t_current < t_final) {
            Scalar dt_step = adjust_step_size(state, dt, t_current);
            state          = static_cast<Derived *>(this)->step(state, dt_step);
            t_current += dt_step;
        }
        Scalar dt_step = t_final - t_current;
        state          = static_cast<Derived *>(this)->step(state, dt_step);
        return state;
    }

    // These functions must be implemented in the derived class
    State step(const State &initial, Scalar dt) {
        return static_cast<Derived *>(this)->step_impl(initial, dt);
    }

protected:
    Function dynamics_function_;

    Scalar adjust_step_size(const State &state, Scalar dt, Scalar current_time) {
        if constexpr (requires { static_cast<Derived *>(this)->adjust_step_size_impl(dt, state, current_time); }) {
            return static_cast<Derived *>(this)->adjust_step_size_impl(dt, state, current_time);
        } else {
            return dt;
        }
    }

    // TODO: Considerations need to be made with this... Ideally we don't want side effects to the global state,
    //       so we will need to provide copies of the subsection of the relevant state and process functions.
    State perform_step(const State &initial, Scalar dt) {
        return static_cast<Derived *>(this)->step_impl(initial, dt);
    }
};

template<typename Derived, typename PositionFunction, typename VelocityFunction>
class BaseSymplecticTraits {
public:
    using pos_Function = PositionFunction;
    using vel_Function = VelocityFunction;

    explicit BaseSymplecticTraits(pos_Function position_function, vel_Function velocity_function)
        : position_function(std::move(position_function)), velocity_function(std::move(velocity_function)) {}

protected:
    pos_Function position_function;
    vel_Function velocity_function;
};

template<typename Derived, typename State, typename Function>
class BaseMultiStepTraits {
public:
    using Derivative = State;
    using Scalar     = State::Scalar;
    using Error      = State;

    // This function must be implemented in the derived class
    virtual State
    step_multistep(const State &initial, Scalar dt, const std::deque<State> &history) = 0;

protected:
    std::deque<State> history_;
};

#endif//BASE_INTEGRATOR_HPP