/**
 * @file state_variant_handler.hpp
 * @author Geoffrey H. Garrett
 * @brief State Variant Handler for ODIN: Optimal Dynamics Integration and Numerics.
 *
 * @details
 * This file provides a mechanism to handle a variant of possible state types
 * in a type-safe manner using std::variant and std::visit.
 *
 * This file is part of ODIN: Optimal Dynamics Integration and Numerics.
 * @copyright 2023 ODIN: Optimal Dynamics Integration and Numerics, All Rights Reserved.
 *
 */

#ifndef STATE_VARIANT_HANDLER_H
#define STATE_VARIANT_HANDLER_H

#include <variant>

/**
 * @brief This class template handles a variant of possible state types.
 *
 * @details
 * A state variant can be any type from the specified States. The StateVariantHandler
 * uses std::visit to dispatch function calls to the correct overload based on the
 * active type within the std::variant at runtime.
 *
 * @tparam States The possible state types that could be in the variant.
 */
template<typename... States>
struct StateVariantHandler {
    using StateVariant = std::variant<States...>;

    /**
      * @brief Applies a function to the state variant in a type-safe manner.
      *
      * @details
      * This function uses std::visit to apply a function to the state variant.
      * The actual function called will depend on the active member of the variant at runtime.
      *
      * @param function The function to apply.
      * @param state The state variant.
      *
      * @return The result of applying the function to the state variant.
      */
    template<typename Function>
    auto visit(Function &&function, const StateVariant &state) const {
        return std::visit(std::forward<Function>(function), state);
    }
};

#endif // STATE_VARIANT_HANDLER_H

/*
 * USAGE EXAMPLES
 *
 * Consider having two state types, `StateA` and `StateB`, and a `StateVariantHandler` for them:
 *
 * struct StateA { void foo() const { std::cout << "StateA's foo\n"; } };
 * struct StateB { void foo() const { std::cout << "StateB's foo\n"; } };
 *
 * using MyHandler = StateVariantHandler<StateA, StateB>;
 *
 * You can create a variant of the states and apply a function to it using the handler:
 *
 * MyHandler::StateVariant myState = StateA{};
 * MyHandler handler;
 *
 * handler.visit([](auto& state) { state.foo(); }, myState);  // Outputs "StateA's foo"
 *
 * If you change the active member of the variant, the function call will adjust accordingly:
 *
 * myState = StateB{};
 * handler.visit([](auto& state) { state.foo(); }, myState);  // Outputs "StateB's foo"
 *
 * The handler thus allows to work with a set of possible state types in a type-safe manner,
 * choosing the correct function overload based on the active state type at runtime.
 */