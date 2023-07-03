#include <iostream>
#include <deque>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <functional>
#include <variant>

#pragma once


template<typename T, typename Container>
struct ContainerPolicy; // forward declaration

/**
 * @brief Policy class for the case when no container is used for history tracking.
 * @tparam T The data type.
 */
template<typename T>
struct ContainerPolicy<T, void> {
    /**
     * @brief Adds data to the container.
     * @param data The data to be added.
     */
    void add(const T &data) {}

    /**
     * @brief Clears the container.
     */
    void clear() {}
};


/**
 * @brief Policy class for sequence containers like std::vector, std::deque, etc.
 * @tparam T The data type.
 * @tparam SeqType The sequence type.
 */
template<typename T, template<typename...> class SeqType>
struct ContainerPolicy<T, SeqType<T>> {
    /**
     * @brief Adds data to the container.
     * @tparam Args Variadic template arguments.
     * @param container The container to be added to.
     * @param args The data to be added.
     */
    template<typename... Args>
    void add(SeqType<T> &container, Args &&...args) {
        static_assert(( std::is_same_v < T, std::remove_reference_t < Args >> && ... ),
        "All arguments must be of the container's value type.");
        ( container.push_back(std::forward<Args>(args)), ... );
    }

    /**
     * @brief Clears the container.
     * @param container The container to be cleared.
     */
    void clear(SeqType<T> &container) {
        container.clear();
    }

    [[nodiscard]] size_t size_in_bytes(const SeqType<T> &container) const {
        // Assuming sizeof gives the complete size for T
        return container.size() * sizeof(T);
    }
};


/**
 * @brief Policy class for associative containers like std::map and std::unordered_map.
 * @tparam Pair The pair type (key, value).
 * @tparam MapType The map type.
 */
template<typename Pair, template<typename...> class MapType>
struct ContainerPolicy<Pair, MapType<typename Pair::first_type, typename Pair::second_type>> {
    using Key = typename Pair::first_type;
    using T = typename Pair::second_type;

    /**
     * @brief Adds data to the container.
     * @param container The container to be added to.
     * @param key The key to be added.
     * @param data The data to be added.
     */
    void add(MapType<Key, T> &container, const Key &key, const T &data) {
        container[key] = data;
    }

    /**
     * @brief Adds a pair to the container.
     * @param container The container to be added to.
     * @param pair The pair to be added.
     */
    void add(MapType<Key, T> &container, const Pair &pair) {
        container[pair.first] = pair.second;
    }

    /**
     * @brief Clears the container.
     * @param container The container to be cleared.
     */
    void clear(MapType<Key, T> &container) {
        container.clear();
    }


    [[nodiscard]] size_t size_in_bytes(const MapType<Key, T> &container) const {
        // Assuming sizeof gives the complete size for key and value types
        return container.size() * (sizeof(Key) + sizeof(T));
    }
};

/**
 * @brief Policy class for array containers working as a sliding window.
 * @tparam T The data type.
 */
template<typename T, std::size_t N>
struct ContainerPolicy<T, std::array<T, N>> {
    std::size_t index_ = 0;

    /**
     * @brief Adds data to the container.
     * @param container The container to be added to.
     * @param data The data to be added.
     */
    void add(std::array<T, N> &container, const T &data) {
        container[index_ % N] = data;
        index_++;
    }

    /**
     * @brief Clears the container.
     * @param container The container to be cleared.
     */
    void clear(std::array<T, N> &container) {
        container.fill(T());
        index_ = 0;
    }

    [[nodiscard]] size_t size_in_bytes(const std::array<T, N> &container) const {
        // Size of the whole array
        return sizeof(container);
    }

    std::array<T, N> get_ordered_container(const std::array<T, N>& container) const {
        std::array<T, N> ordered_container;
        for (std::size_t i = 0; i < N; ++i) {
            ordered_container[i] = container[(i + index_) % N];
        }
        return ordered_container;
    }
};