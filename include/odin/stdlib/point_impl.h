#ifndef ODIN_STDLIB_POINT_IMPL_H
#define ODIN_STDLIB_POINT_IMPL_H

#include <odin/core/point.h>

namespace odin::shape {


//    // Specialization for std::array<T, N>
//    template<typename T, std::size_t N>
//    struct point_series_type_trait<std::array<T, N>> {
//        using type = std::array<std::array<T, N>, N>; // Two-dimensional array
//    };


// Array squared_distance
template<Point T, std::size_t N>
    requires std::same_as<T, std::array<scalar_type<T>, N>>
constexpr scalar_type<T> squared_distance(const T &a, const T &b) {
#ifdef __cpp_lib_transform_reduce
    return std::transform_reduce(
            a.begin(), a.end(), b.begin(), scalar_type<T>{0}, std::plus{},
            [](auto a, auto b) { return std::pow(a - b, 2); });
#else
    scalar_type<T> sum = 0;
    for (std::size_t i = 0; i < N; ++i) { sum += std::pow(a[i] - b[i], 2); }
    return sum;
#endif
}

// Array squared_norm
template<Point T, std::size_t N>
    requires std::same_as<T, std::array<scalar_type<T>, N>>
constexpr scalar_type<T> squared_norm(const T &a) {
#ifdef __cpp_lib_transform_reduce
    return std::transform_reduce(a.begin(), a.end(), scalar_type<T>{0},
                                 std::plus{},
                                 [](auto a) { return std::pow(a, 2); });
#else
    scalar_type<T> sum = 0;
    for (std::size_t i = 0; i < N; ++i) { sum += std::pow(a[i], 2); }
    return sum;
#endif
}

// Array element_wise_division
template<Point T, std::size_t N>
    requires std::same_as<T, std::array<scalar_type<T>, N>>
constexpr T element_wise_division(const T &a, const T &b) {
    T result;
    std::ranges::transform(a, b, result.begin(), std::divides{});
    return result;
}

template<Point T, std::size_t N>
    requires std::same_as<T, std::array<scalar_type<T>, N>>
constexpr scalar_type<T> dot_product(const T &a, const T &b) {
#ifdef __cpp_lib_transform_reduce
    return std::transform_reduce(a.begin(), a.end(), b.begin(),
                                 scalar_type<T>{0}, std::plus{},
                                 std::multiplies{});
#else
    scalar_type<T> sum = 0;
    for (std::size_t i = 0; i < N; ++i) { sum += a[i] * b[i]; }
    return sum;
#endif
}

// Array zero
template<Point T, std::size_t N>
    requires std::same_as<T, std::array<scalar_type<T>, N>>
constexpr T zero() {
    return T{};
}

// Array Euclidean norm
template<typename T, std::size_t N>
    requires std::same_as<T, std::array<scalar_type<T>, N>>
constexpr scalar_type<T> normalize(const T &a) {
    return std::sqrt(std::accumulate(a.begin(), a.end(),
                                     static_cast<scalar_type<T>>(0),
                                     [](const auto &sum, const auto &element) {
                                         return sum + element * element;
                                     }));
}
}// namespace odin::traits

#endif// ODIN_STDLIB_POINT_IMPL_H