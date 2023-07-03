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

template<typename T>
struct ContainerPolicy<T, void> {
    void add(const T &data) {}

    void clear() {}
};

template<typename T, template<typename...> class SeqType>
struct ContainerPolicy<T, SeqType<T>> {
    template<typename... Args>
    void add(SeqType<T> &container, Args &&...args) {
        static_assert(( std::is_same_v < T, std::remove_reference_t < Args >> && ... ),
        "All arguments must be of the container's value type.");
        ( container.push_back(std::forward<Args>(args)), ... );
    }

    void clear(SeqType<T> &container) {
        container.clear();
    }

    [[nodiscard]] size_t size_in_bytes(const SeqType<T> &container) const {
        // Assuming sizeof gives the complete size for T
        return container.size() * sizeof(T);
    }
};

template<typename Pair, template<typename...> class MapType>
struct ContainerPolicy<Pair, MapType<typename Pair::first_type, typename Pair::second_type>> {
    using Key = typename Pair::first_type;
    using T = typename Pair::second_type;

    void add(MapType<Key, T> &container, const Key &key, const T &data) {
        container[key] = data;
    }

    void add(MapType<Key, T> &container, const Pair &pair) {
        container[pair.first] = pair.second;
    }

    void clear(MapType<Key, T> &container) {
        container.clear();
    }

    [[nodiscard]] size_t size_in_bytes(const MapType<Key, T> &container) const {
        return container.size() * (sizeof(Key) + sizeof(T));
    }
};

template<typename T, std::size_t N>
struct ContainerPolicy<T, std::array<T, N>> {
    std::size_t index_ = 0;

    void add(std::array<T, N> &container, const T &data) {
        container[index_ % N] = data;
        index_++;
    }

    void clear(std::array<T, N> &container) {
        container.fill(T());
        index_ = 0;
    }

    [[nodiscard]] size_t size_in_bytes(const std::array<T, N> &container) const {
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