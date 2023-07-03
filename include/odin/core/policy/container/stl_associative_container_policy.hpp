#include <set>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include "container_policy.hpp"

// Create type traits to distinguish between map-like and set-like containers
template<typename T>
struct is_map_container : std::false_type {
};

template<typename Key, typename T, typename... Args>
struct is_map_container<std::map<Key, T, Args...>> : std::true_type {
};

template<typename Key, typename T, typename... Args>
struct is_map_container<std::multimap<Key, T, Args...>> : std::true_type {
};

template<typename Key, typename T, typename... Args>
struct is_map_container<std::unordered_map<Key, T, Args...>> : std::true_type {
};

template<typename Key, typename T, typename... Args>
struct is_map_container<std::unordered_multimap<Key, T, Args...>> : std::true_type {
};

template<typename T>
inline constexpr bool is_map_container_v = is_map_container<T>::value;


template<typename T>
struct is_set_container : std::false_type {
};

template<typename Key, typename... Args>
struct is_set_container<std::set<Key, Args...>> : std::true_type {
};

template<typename Key, typename... Args>
struct is_set_container<std::multiset<Key, Args...>> : std::true_type {
};

template<typename Key, typename... Args>
struct is_set_container<std::unordered_set<Key, Args...>> : std::true_type {
};

template<typename Key, typename... Args>
struct is_set_container<std::unordered_multiset<Key, Args...>> : std::true_type {
};

template<typename T>
inline constexpr bool is_set_container_v = is_set_container<T>::value;

// Create type traits to distinguish between associative and unordered associative containers
template<typename T>
struct is_associative_container : std::false_type {
};

template<typename Key, typename T, typename... Args>
struct is_associative_container<std::map<Key, T, Args...>> : std::true_type {
};

template<typename Key, typename T, typename... Args>
struct is_associative_container<std::multimap<Key, T, Args...>> : std::true_type {
};

template<typename Key, typename... Args>
struct is_associative_container<std::set<Key, Args...>> : std::true_type {
};

template<typename Key, typename... Args>
struct is_associative_container<std::multiset<Key, Args...>> : std::true_type {
};

template<typename T>
inline constexpr bool is_associative_container_v = is_associative_container<T>::value;

template<typename T>
struct is_unordered_associative_container : std::false_type {
};

template<typename Key, typename T, typename... Args>
struct is_unordered_associative_container<std::unordered_map<Key, T, Args...>> : std::true_type {
};

template<typename Key, typename T, typename... Args>
struct is_unordered_associative_container<std::unordered_multimap<Key, T, Args...>> : std::true_type {
};

template<typename Key, typename... Args>
struct is_unordered_associative_container<std::unordered_set<Key, Args...>> : std::true_type {
};

template<typename Key, typename... Args>
struct is_unordered_associative_container<std::unordered_multiset<Key, Args...>> : std::true_type {
};

template<typename T>
inline constexpr bool is_unordered_associative_container_v = is_unordered_associative_container<T>::value;

template<typename T>
struct is_pair : std::false_type {};

template<typename T, typename U>
struct is_pair<std::pair<T, U>> : std::true_type {};

template<typename T>
inline constexpr bool is_pair_v = is_pair<T>::value;

#include "container_policy.hpp"


/**
 * @brief Computes the dynamic size of an object.
 *
 * This function is used to compute the dynamic size of an object. It is specialized to handle
 * several common types of STL containers including std::string, std::vector, and std::deque.
 * For containers, it computes the total size of all contained objects. It also handles nested
 * containers by recursive calls to `dynamic_size`, up to a specified depth.
 *
 * @param object The object whose dynamic size is to be computed.
 * @param depth The maximum recursion depth. By default, this is 0 (i.e., recursion is disabled).
 * @return The total dynamic size in bytes of the object.
 *
 * @note This function introduces overhead due to recursion for nested containers.
 * It may not accurately measure the memory usage for complex data structures due to potential padding.
 * For types not handled by this function, it returns 0.
 */
template<typename T>
size_t dynamic_size(const T &object, int depth = 0) {
    if constexpr (std::is_same_v<T, std::string>) {
        // For std::string, return the number of characters.
        return object.size();
    } else if constexpr (std::is_same_v<T, std::vector<typename T::value_type>>) {
        // For std::vector, return the total size of all elements.
        size_t total_size = 0;
        if (depth > 0) {
            for (const auto &element: object) {
                total_size += dynamic_size(element, depth - 1);
            }
        } else {
            total_size = object.size() * sizeof(typename T::value_type);
        }
        return total_size;
    } else if constexpr (std::is_same_v<T, std::deque<typename T::value_type>>) {
        // For std::deque, return the total size of all elements.
        size_t total_size = 0;
        if (depth > 0) {
            for (const auto &element: object) {
                total_size += dynamic_size(element, depth - 1);
            }
        } else {
            total_size = object.size() * sizeof(typename T::value_type);
        }
        return total_size;
    } else {
        // For other types, return 0.
        return 0;
    }
}

#include <cereal/types/polymorphic.hpp>


/**
 * @class ContainerPolicy
 *
 * @brief Policy class for handling map containers.
 *
 * This class is specialized for map containers, providing functionality that is
 * specifically tuned to the characteristics of map containers.
 *
 * @tparam MapType The type of the map container.
 * @tparam Key The type of the key used in the map.
 * @tparam T The type of the value stored in the map.
 * @tparam Args Additional template arguments.
 *
 * @note This class is currently written using C++20's "requires" clause for
 *       simplicity and improved error messages. However, performance should be
 *       identical to the equivalent C++17 version that uses std::enable_if_t.
 *       If needed, you can switch to the C++17 version by replacing the "requires"
 *       clause with the commented out code:
 *
 *       @code{.cpp}
 *       template<
 *           template<typename, typename, typename...> class MapType,
 *           typename Key, typename T, typename... Args,
 *           typename = std::enable_if_t<is_map_container_v<MapType<Key, T, Args...>>>
 *       >
 *       struct ContainerPolicy<MapType<Key, T, Args...>>
 *       @endcode
 */
template<template<typename, typename, typename...> class MapType, typename Key, typename T, typename... Args> requires is_map_container_v<MapType<Key, T, Args...>>
struct ContainerPolicy<MapType<Key, T, Args...>>
        : public ContainerPolicyBase<ContainerPolicy<MapType<Key, T, Args...>>, MapType<Key, T, Args...>> {
    using Container = MapType<Key, T, Args...>;
    using Pair = std::pair<const Key, T>;

    // Adding a pair
    template<typename Arg>
    std::enable_if_t<is_pair_v<std::decay_t<Arg>>, void>
    add_impl(Container &container, Arg &&arg) {
        container.insert(std::forward<Arg>(arg));
    }

    // Adding key and value separately
    template<typename K, typename V>
    std::enable_if_t<!std::is_same_v<std::decay_t<K>, Pair> && std::is_same_v<std::decay_t<V>, T>, void>
    add_impl(Container &container, K &&key, V &&value) {
        container.insert(std::make_pair(std::forward<K>(key), std::forward<V>(value)));
    }

    void clear_impl(Container &container) {
        container.clear();
    }

    template<typename Archive>
    void save_impl(Archive & ar, const Container &container) const {
        ar(cereal::make_size_tag(static_cast<cereal::size_type>(container.size()))); // number of elements
        for(auto& [key, value]: container) {
            ar(cereal::make_map_item(key, value));
        }
    }

    template<typename Archive>
    void load_impl(Archive& ar, Container& container) {
        cereal::size_type size;
        ar(cereal::make_size_tag(size));

        container.clear();
        for(cereal::size_type i = 0; i < size; ++i) {
            Key key;
            T value;
            ar(key, value);
            container.insert({key, value});
        }
    }
    /**
      * @brief Computes the total size in bytes of the container.
      *
      * This function computes the total size of the container in bytes. It includes the size of each key-value pair
      * and the overhead of each node in a tree-based map or similar container. It can optionally calculate the
      * dynamic size of elements if they have the dynamic_size function available.
      *
      * Each node in a tree-based container (such as std::map) typically has an overhead for:
      * - The key-value pair (sizeof(Pair))
      * - Two child pointers (2 * sizeof(void*))
      * - The parent pointer (sizeof(void*))
      * - The color property (sizeof(bool))
      *
      * @param container Reference to the container whose size is to be computed.
      * @return The total size in bytes of the container.
      *
      * @note This function can introduce performance overhead if the dynamic_size function is used.
      * It does not handle potential padding inside complex data structures.
      * It assumes that all dynamic memory used by elements can be measured with dynamic_size.
    */
    [[nodiscard]] size_t size_in_bytes_impl(const Container &container) const {
        size_t dynamic_key_size = 0;
        size_t dynamic_value_size = 0;

        // Compute the dynamic size of keys, if available.
        if constexpr (std::is_invocable_v<decltype(dynamic_size<Key>), const Key &>) {
            for (const auto &[key, _]: container) {
                dynamic_key_size += dynamic_size(key);
            }
        }

        // Compute the dynamic size of values, if available.
        if constexpr (std::is_invocable_v<decltype(dynamic_size<T>), const T &>) {
            for (const auto &[_, value]: container) {
                dynamic_value_size += dynamic_size(value);
            }
        }

        // Compute the total size, including the size of each element and the overhead of each node.
        return container.size() * (sizeof(Pair) + 3 * sizeof(void *) + sizeof(bool)) + dynamic_key_size +
               dynamic_value_size;
    }

};

/**
 * @class ContainerPolicy
 *
 * @brief Policy class for handling set containers.
 *
 * This class is specialized for set containers, providing functionality that is
 * specifically tuned to the characteristics of set containers.
 *
 * @tparam SetType The type of the set container.
 * @tparam Key The type of the key used in the set.
 * @tparam Args Additional template arguments.
 *
 * @note This class is currently written using C++20's "requires" clause for
 *       simplicity and improved error messages. However, performance should be
 *       identical to the equivalent C++17 version that uses std::enable_if_t.
 *       If needed, you can switch to the C++17 version by replacing the "requires"
 *       clause with the commented out code:
 *
 *       @code{.cpp}
 *       template<
 *           template<typename, typename...> class SetType,
 *           typename Key, typename... Args,
 *           typename = std::enable_if_t<is_set_container_v<SetType<Key, Args...>>>
 *       >
 *       struct ContainerPolicy<SetType<Key, Args...>>
 *       @endcode
 */
template<template<typename, typename, typename...> class SetType, typename Key, typename... Args> requires is_set_container_v<SetType<Key, Args...>>
struct ContainerPolicy<SetType<Key, Args...>>
        : public ContainerPolicyBase<ContainerPolicy<SetType<Key, Args...>>, SetType<Key, Args...>> {
    using Container = SetType<Key, Args...>;

    // Adding a key
    template<typename K>
    void add_impl(Container &container, K &&key) {
        container.insert(std::forward<K>(key));
    }

    void clear_impl(Container &container) {
        container.clear();
    }

    /**
     * @brief Computes the total size in bytes of the container.
     *
     * This function computes the total size of the container in bytes. It includes the size of each key
     * and the overhead of each node in a tree-based set or similar container. It can optionally calculate the
     * dynamic size of keys if they have the dynamic_size function available.
     *
     * Each node in a tree-based container (such as std::set) typically has an overhead for:
     * - The key (sizeof(Key))
     * - Two child pointers (2 * sizeof(void*))
     * - The parent pointer (sizeof(void*))
     * - The color property (sizeof(bool))
     *
     * @param container Reference to the container whose size is to be computed.
     * @return The total size in bytes of the container.
     *
     * @note This function can introduce performance overhead if the dynamic_size function is used.
     * It does not handle potential padding inside complex data structures.
     * It assumes that all dynamic memory used by keys can be measured with dynamic_size.
     */
    [[nodiscard]] size_t size_in_bytes_impl(const Container &container) const {
        size_t dynamic_key_size = 0;

        // Compute the dynamic size of keys, if available.
        if constexpr (std::is_invocable_v<decltype(dynamic_size<Key>), const Key &>) {
            for (const auto &key: container) {
                dynamic_key_size += dynamic_size(key);
            }
        }

        // Compute the total size, including the size of each element and the overhead of each node.
        return container.size() * (sizeof(Key) + 3 * sizeof(void *) + sizeof(bool)) + dynamic_key_size;
    }

};

