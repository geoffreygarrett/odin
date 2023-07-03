
// x{i+1} = f(x{i}, t{i}) + q{i}, only if noise is considered

#include <tuple>

/**
 * @brief Template class for handling state objects by splitting and combining them.
 *
 * @tparam State The type of the primary state object.
 * @tparam States List of additional state types to be handled.
 */
template<typename State, typename... States>
class StateHandler {
public:
    /**
     * @brief Splits a state object into multiple tuples by invoking the split() member function of each state type.
     *
     * @tparam T The type of the state object.
     * @param state The state object to be split.
     * @return A tuple containing the split state components.
     */
    template<typename T>
    static auto split(const T &state) {
        return std::tuple_cat(State::split(state), StateHandler<States...>::split(state));
    }

    /**
     * @brief Combines a state object with multiple tuples by invoking the combine() member function of each state type.
     *
     * @tparam T The type of the state object.
     * @tparam TupleTypes The types of the tuples to be combined.
     * @param state The state object to be combined.
     * @param tuples The tuples to be combined with the state object.
     */
    template<typename T, typename... TupleTypes>
    static void combine(T &state, TupleTypes &&... tuples) {
        (State::combine(state, std::forward<TupleTypes>(tuples)), ...);
        StateHandler<States...>::combine(state, std::forward<TupleTypes>(tuples)...);
    }
};
