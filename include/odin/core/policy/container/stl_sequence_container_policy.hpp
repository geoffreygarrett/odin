#include "container_policy.hpp"

/**
 * @brief Generic ContainerPolicy for STL sequence containers.
 *
 * This policy class defines common operations for STL sequence containers like std::vector, std::list, etc.
 * It includes methods for adding elements, clearing the container, and calculating the total size of the
 * container in bytes. The size calculation takes into account both the static size of elements and, optionally,
 * their dynamic size.
 *
 * @tparam ContainerType The type of the container.
 * @tparam T The type of the elements in the container.
 * @tparam Args Any additional template arguments for the container type.
 */
template<template<typename, typename...> class ContainerType, typename T, typename... Args>
struct ContainerPolicy<ContainerType<T, Args...>>
        : public ContainerPolicyBase<ContainerPolicy<ContainerType<T, Args...>>, ContainerType<T, Args...>> {
    using Container = ContainerType<T, Args...>;

    // Add one or more elements to the container.
    template<typename... AddArgs>
    void add_impl(Container &container, AddArgs &&... args) {
        ( container.push_back(std::forward<AddArgs>(args)), ... );
    }

    // Clear the container.
    void clear_impl(Container &container) {
        container.clear();
    }

    /**
      * @brief Computes the total size in bytes of the container.
      *
      * This function computes the total size of the container in bytes, including both the static size of
      * elements and their dynamic size, if they have the dynamic_size function available.
      *
      * @param container The container whose size is to be computed.
      * @return The total size in bytes of the container.
      *
      * @note This function introduces overhead due to the calculation of dynamic size.
      * It may not accurately measure the memory usage for complex data structures due to potential padding.
      * For elements not handled by dynamic_size, it only includes their static size.
    */
    [[nodiscard]] size_t size_in_bytes_impl(const Container &container) const {
        size_t total_size = 0;
        if constexpr (std::is_invocable_v<decltype(dynamic_size<T>), const T &>) {
            for (const auto &element: container) {
                total_size += dynamic_size(element);
            }
        } else {
            total_size = container.size() * sizeof(T);
        }
        return total_size;
    }
};

/**
 * @brief ContainerPolicy for std::array.
 *
 * This policy class defines common operations for std::array. It includes methods for adding elements,
 * clearing the array, calculating the total size of the array in bytes, and obtaining an ordered version of
 * the array. The size calculation takes into account both the static size of elements and their dynamic size,
 * if they have the dynamic_size function available.
 *
 * @tparam T The type of the elements in the array.
 * @tparam N The size of the array.
 */
template<typename T, std::size_t N>
struct ContainerPolicy<std::array<T, N>> {
    std::size_t index_ = 0;

    // Add an element to the array.
    void add(std::array<T, N> &container, const T &data) {
        container[index_ % N] = data;
        index_++;
    }

    // Clear the array.
    void clear(std::array<T, N> &container) {
        container.fill(T());
        index_ = 0;
    }

    /**
      * @brief Computes the total size in bytes of the array.
      *
      * This function computes the total size of the array in bytes, including both the static size of
      * elements and their dynamic size, if they have the dynamic_size function available.
      *
      * @param container The array whose size is to be computed.
      * @return The total size in bytes of the array.
      *
      * @note This function introduces overhead due to the calculation of dynamic size.
      * It may not accurately measure the memory usage for complex data structures due to potential padding.
      * For elements not handled by dynamic_size, it only includes their static size.
    */
    [[nodiscard]] size_t size_in_bytes(const std::array<T, N> &container) const {
        size_t total_size = 0;
        if constexpr (std::is_invocable_v<decltype(dynamic_size<T>), const T &>) {
            for (const auto &element: container) {
                total_size += dynamic_size(element);
            }
        } else {
            total_size = sizeof(container);
        }
        return total_size;
    }

    // Get an ordered version of the array.
    std::array<T, N> get_ordered_container(const std::array<T, N> &container) const {
        std::array<T, N> ordered_container;
        for (std::size_t i = 0; i < N; ++i) {
            ordered_container[i] = container[(i + index_) % N];
        }
        return ordered_container;
    }
};