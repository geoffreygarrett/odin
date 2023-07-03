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


/**
 * @class StateHandler
 *
 * @brief A handler for managing the states of celestial bodies.
 *
 * This class is a generic handler for managing the states of various celestial bodies.
 * It's designed to work with a variety of state types, provided as template parameters.
 *
 * @tparam States Types of the celestial body states.
 *
 * @note This class uses several modern C++ features for more efficient and expressive code:
 *       - Parameter packs (from C++11): for handling an arbitrary number of state types.
 *       - Tuple and apply function (from C++17): for storing states and operating on them.
 *       - Perfect forwarding (from C++11) in combine() and recombine(): to avoid unnecessary copies and
 *         improve performance. These methods expect rvalue references and will consume their arguments.
 *         The caller should ensure that any necessary data is saved prior to calling these functions.
 *       - Auto type deduction (from C++14) in get(): to avoid writing out the complicated return type.
 *       - Structured bindings (from C++17) in recombine(): for cleaner unpacking of the state tuple.
 */
template<StateConcept... States>
class StateHandler {
public:
    // Tuple that stores individual states. Each type in the States parameter pack represents a different state, and
    // this tuple contains an instance of each. These are state traits or properties and are not owned by the tuple.
    // Instead, they offer a non-owning view into the actual state data which is stored contiguously elsewhere.
    // Note that each state type (in States...) must have an associated Eigen type to represent its data.
    // WARNING: Care must be taken to manage the lifetime of the actual state data.
    // The data must outlive the views provided here, and should not be relocated.
    // Failure to ensure this may lead to undefined behavior.
    using state_tuple_t = std::tuple<States...>;

    // A static constexpr variable to calculate the total dimension of the state vector.
    // It uses a fold expression (a feature from C++17) over parameter pack of States.
    // The fold expression (... + States::state_dimension) sums up the 'state_dimension' of each State in the pack.
    // Each State should have a static constexpr member 'state_dimension' specifying its dimension.
    // The final result 'state_dimension' gives the total dimension of the state vector represented by the CompositeState.
    static constexpr std::size_t state_dimension = (... + States::state_dimension);

    // now define the scalar type using first state in the tuple
    using scalar_t = typename std::tuple_element<0, state_tuple_t>::type::scalar_type;

    // Retrieve all State::value matrices, and float/double and combine them into a single Eigen vector at compile
    // time. This is used to create a single Eigen vector from the state_tuple.
    // This is a helper function for the constructor.
    // The function uses a fold expression (a feature from C++17) over parameter pack of States.
    using state_vector_t = Eigen::Matrix<scalar_t, state_dimension, 1>;

    // static assert to ensure that the state_dimension is not zero.
    static_assert(state_dimension != 0, "State dimension cannot be zero.");

    // static assert to ensure all scalar values are the same in the tuple
    // (i.e., all states use the same scalar type).
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

// NOTE: The DynamicStateHandler is highly experimental and not yet ready for use. This is more of a proof of concept in
// case we want to implement dynamic states in the future.
// NOTE: This implementation and the following one, both assume the existence of "Entity" which the static part does not.
// I find this preferable, but it may be difficult to stay agnostic to the existence of Entity, if we add dynamic states.
// NOTE: If we started using dynamic states, we would need to consider memory management, as we would need to allocate
// memory for the states at runtime. This is not a problem for the static states, as they are allocated at compile time.
#ifdef ENABLE_DYNAMIC_STATES
    // Dynamic State Handler
    class DynamicStateHandler {
        std::unordered_map<std::string, std::shared_ptr<Entity>> dynamic_entities;

    public:
        template <typename T, typename... Args>
        void emplace(Args&&... args) {
            dynamic_entities.emplace(T::name, std::make_shared<T>(std::forward<Args>(args)...));
        }

        Entity& get(const std::string& name) {
            return *dynamic_entities.at(name);
        }
    } dynamic;

    // Overload of get function for runtime entities.
    Entity& get(const std::string& name) {
        return dynamic.get(name);
    }
#endif
// NOTE: This is the alternative implementation of dynamic states, it may be a little less distinct for users as to
// when we are using static or dynamic states, but it *might* be more intuitive to use.
#ifdef ENABLE_DYNAMIC_STATES_ALTERNATIVE
    /**
     * @brief Adds a new runtime entity to the system.
     *
     * @param id The ID of the new entity.
     * @param states The state tuple of the new entity.
     */
    static void add_entity(const std::string& id, state_tuple&& states) {
        runtime_entities.emplace(id, std::move(states));
    }

    /**
     * @brief Removes a runtime entity from the system.
     *
     * @param id The ID of the entity to remove.
     */
    static void remove_entity(const std::string& id) {
        runtime_entities.erase(id);
    }

    /**
     * @brief Gets a state of a runtime entity by its type.
     *
     * @tparam T The type of the state to get.
     * @param id The ID of the entity.
     * @return A reference to the state of type T.
     */
    template<typename T>
    [[nodiscard]] static auto &get_runtime(const std::string& id) {
        return std::get<T>(runtime_entities[id]);
    }
#endif


};

#endif //STATE_HANDLER_HPP