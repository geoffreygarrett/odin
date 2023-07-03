#ifndef STATE_HANDLER_HPP
#define STATE_HANDLER_HPP

#include <iostream>
#include <tuple>
#include <type_traits>
#include <Eigen/Core>
#include "state.hpp"


template<typename... States>
struct StateSize {
    static constexpr int value = (... + States::state_dimension);
};



template<StateConcept... States>
class StateHandler {
public:

    using state_tuple_t = std::tuple<States...>;
    static constexpr std::size_t state_dimension = (... + States::state_dimension);
    using scalar_t = typename std::tuple_element<0, state_tuple_t>::type::scalar_type;
    using state_vector_t = Eigen::Matrix<scalar_t, state_dimension, 1>;
    static_assert(state_dimension != 0, "State dimension cannot be zero.");
    static_assert(
            (std::is_same_v<typename States::scalar_type, typename States::scalar_type> && ...),
            "All states must use the same scalar type. (not STRICTLY necessary, but rather be consistent am I right?)");

//private:
     state_tuple_t state_tuple;
     state_vector_t state_vector;

    // Initialize the state views from the provided state data and initialize the global state vector
//    explicit StateHandler(States... states) : state_tuple(std::move(states)...) {
//        init_global_state(); // initialize the global state vector from the provided state data, in the tuples
//        init_state_views(); // initialize the state views to point to the correct parts of the global state vector
//    }

    // Initialize the global state vector from the provided state data
    void init_global_state() {
        std::size_t idx = 0;
        // Use a lambda expression to iterate over the states in the tuple
        auto init_state = [&](auto &state) {
            // Extract the data from the state and assign it to the global state vector
            for (int i = 0; i < state.state_dimension; ++i) {
                state_vector(idx + i) = state.value(i);
            }
            // Advance the index
            idx += state.state_dimension;
        };
        std::apply([&](auto &... states) { (init_state(states), ...); }, state_tuple);
    }

    // Initialize the state views to point to the correct parts of the global state vector
    void init_state_views() {
        std::size_t idx = 0;
        // Use a lambda expression to iterate over the states in the tuple
        auto init_view = [&](auto &state) {
            // Replace the state's data with a view into the global state vector
            state.value = state_vector.segment(idx, state.state_dimension);
            // Advance the index
            idx += state.state_dimension;
        };
        std::apply([&](auto &... states) { (init_view(states), ...); }, state_tuple);
    }

    template<typename T>
    [[nodiscard]] static auto &get(state_tuple_t &states) {
        return std::get<T>(states);
    }

    template<typename T>
    [[nodiscard]] static const auto &get(const state_tuple_t &states) {
        return std::get<T>(states);
    }

    template<typename ID, typename Tuple, std::size_t... I>
    auto get_states_of_entity_impl(Tuple &&t, std::index_sequence<I...>) {
        return std::tuple_cat(
                (std::is_same_v<ID, typename std::tuple_element_t<I, std::decay_t<Tuple>>::id_type> ?
                 std::make_tuple(std::get<I>(t)) : std::tuple<>())...);
    }

    template<typename ID, typename Tuple>
    auto get_states_of_entity(Tuple &&t) {
        return get_states_of_entity_impl<ID>(
                std::forward<Tuple>(t),
                std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{});
    }

    [[nodiscard]] static auto split(state_tuple_t &states) {
        return std::apply([](auto &... state) {
            return std::make_tuple(std::addressof(state)...);
        }, states);
    }

    static void combine(state_tuple_t &old_states, state_tuple_t &&new_states) {
        old_states = std::move(new_states);
    }

    template<typename... SpecificStates>
    static void recombine(state_tuple_t &old_states, std::tuple<SpecificStates...> &&new_states) {
        std::apply([&old_states](auto &&... specific_state) {
            ((std::get<std::remove_reference_t<decltype(specific_state)>>(
                    old_states) = std::forward<decltype(specific_state)>(specific_state)), ...);
        }, std::move(new_states));
    }

    template<int DerivativeOrder, typename Tuple, std::size_t... I>
    auto get_states_of_order_impl(Tuple &&t, std::index_sequence<I...>) {
        return std::tuple_cat(
                (std::tuple_element_t<I, std::decay_t<Tuple>>::derivative_order == DerivativeOrder ?
                 std::make_tuple(std::get<I>(t)) : std::tuple<>())...);
    }

    template<int DerivativeOrder, typename Tuple>
    auto get_states_of_order(Tuple &&t) {
        return get_states_of_order_impl<DerivativeOrder>(
                std::forward<Tuple>(t),
                std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple>>>{});
    }

};

#endif //STATE_HANDLER_HPP