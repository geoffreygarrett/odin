/**
 * @file
 * @brief Header file defining the Entity system.
 */

// Old minimal implementation
////struct Entity {
////    static constexpr bool is_id = true;
////};
////
////#define DEFINE_ENTITY(NAME) struct NAME : public Entity<NAME> {\
////    static constexpr const char* name = #NAME;\
////}
////
////// The bodies are identified by their types.
////struct EARTH : public Entity {
////    static constexpr const char* name = "EARTH";
////};
////
////struct DELFI : public Entity {
////    static constexpr const char* name = "DELFI";
////};

#include <string>
#include <type_traits>
#include <utility>

// --- Entity Base Classes ---

/**
 * @class StaticEntityTag
 * @brief Tag class for static entities.
 */
struct StaticEntityTag {
};

/**
 * @class RuntimeEntityTag
 * @brief Tag class for runtime entities.
 */
struct RuntimeEntityTag {
};

/**
 * @class Entity
 * @brief Base class for static entities.
 *
 * This class is intended to be used through CRTP, with the derived class as the template parameter.
 * It provides the static constexpr member variables 'is_id' and 'is_runtime', as well as the alias 'EntityTag' for StaticEntityTag.
 * It also provides a constexpr function 'name()' to access the name of the entity, which is defined in the derived class.
 */
template<typename Derived>
struct Entity {
    static constexpr bool is_id = true;
    static constexpr bool is_runtime = false;
    using EntityTag = StaticEntityTag;

    constexpr auto name() const {
        return Derived::name;
    }
};

/**
 * @brief Macro for defining static entity types.
 *
 * This macro defines a struct with the given name that derives from Entity using CRTP,
 * and provides a static constexpr member variable 'name' with the same value as the name of the struct.
 */
#define DEFINE_ENTITY(NAME) struct NAME : public Entity<NAME> {\
    static constexpr const char* name = #NAME;\
}

/**
 * @class RuntimeEntity
 * @brief Class for runtime entities.
 *
 * This class provides the static constexpr member variables 'is_id' and 'is_runtime', as well as the alias 'EntityTag' for RuntimeEntityTag.
 * It also provides a non-constexpr function 'name()' to access the name of the entity, which can be changed at runtime.
 */
struct RuntimeEntity {
    static constexpr bool is_id = true;
    static constexpr bool is_runtime = true;
    using EntityTag = RuntimeEntityTag;
    std::string dynamic_name;

    explicit RuntimeEntity(std::string name) : dynamic_name(std::move(name)) {}

    [[nodiscard]] std::string name() const { return dynamic_name; }
};

// --- Usage Example ---

/**
 * @brief Usage example of the Entity system.
 *
 * This function creates instances of the static entity types and a runtime entity,
 * and prints out the name of each one.
 */
/*

#include <state/entity.hpp>

DEFINE_ENTITY(STATIC_ENTITY);

int main() {
    // Static entity usage
    STATIC_ENTITY static_entity;
    std::cout << "Static Entity: " << static_entity.name() << std::endl; // prints "STATIC_ENTITY"

    // Runtime entity usage
    RUNTIME_ENTITY runtime_entity("RUNTIME_ENTITY");
    std::cout << "Initial Runtime Entity: " << runtime_entity.name() << std::endl; // prints "RUNTIME_ENTITY"

    // Changing the name at runtime
    runtime_entity.dynamic_name = "NEW_RUNTIME_ENTITY_NAME";
    std::cout << "Modified Runtime Entity: " << runtime_entity.name() << std::endl; // prints "NEW_RUNTIME_ENTITY_NAME"

    return 0;
}
*/