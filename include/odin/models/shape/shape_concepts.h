#ifndef ODIN_SHAPE_CONCEPTS_H
#define ODIN_SHAPE_CONCEPTS_H

#include <Eigen/Core>
//#include <unsupported/Eigen/CXX11/Tensor>
#include "unsupported/Eigen/CXX11/Tensor"
#include <array>
#include <bit>
#include <cmath>
#include <numeric>
#include <ranges>
#include <type_traits>

////--- Temporary home of Ray ---//
//template<Point P>
//class Ray {
//public:
//    using point_type = P;
//
//    Ray(const point_type &source, const point_type &direction)
//        : m_source(source), m_direction(direction) {}
//
//    point_type get_source() const { return m_source; }
//    point_type get_direction() const { return m_direction; }
//
//private:
//    point_type m_source;
//    point_type m_direction;
//};


namespace odin::shape {


// Scalar type trait
template<typename T>
struct scalar_type_trait {
    using type = typename T::Scalar;
};

template<typename T, std::size_t N>
struct scalar_type_trait<std::array<T, N>> {
    using type = T;
};

template<typename T>
using scalar_type = typename scalar_type_trait<T>::type;

// Point series type trait
template<typename T>
struct point_series_type_trait {
    using type = std::vector<typename T::point_type>;
};

template<typename T, int N>
struct point_series_type_trait<Eigen::Matrix<T, N, 1>> {
    using type = Eigen::Matrix<T, N, Eigen::Dynamic>;
};

template<typename T>
using point_series_type = typename point_series_type_trait<T>::type;


//template<typename T>
//struct point_grid_type_trait<Eigen::Matrix<T, 2, 1>, 3> {
//    // For a 3D grid of 2D Eigen points, use a 3D Eigen::Tensor
//    // The resulting type is a 3D tensor of 2D Eigen::Vector values
//    using type = Eigen::Tensor<Eigen::Matrix<T, 2, 1>, 3>;
//};
//
//template<typename T>
//struct point_grid_type_trait<Eigen::Matrix<T, 3, 1>, 3> {
//    // For a 3D grid of 3D Eigen points, use a 3D Eigen::Tensor
//    // The resulting type is a 3D tensor of 3D Eigen::Vector values
//    using type = Eigen::Tensor<Eigen::Matrix<T, 3, 1>, 3>;
//};


// Function prototypes for any Point
template<typename T>
constexpr scalar_type<T> squared_distance(const T &a, const T &b);

template<typename T>
constexpr scalar_type<T> squared_norm(const T &a);

template<typename T>
constexpr T element_wise_division(const T &a, const T &b);

template<typename T>
constexpr T zero();

template<typename T>
constexpr scalar_type<T> dot_product(const T &a, const T &b);

template<typename T>
constexpr T normalize(const T &a);

// Concept for Point
template<typename T>
concept Point = requires(T a, T b) {
    { squared_distance(a, b) } -> std::same_as<scalar_type<T>>;
    { squared_norm(a) } -> std::same_as<scalar_type<T>>;
    { element_wise_division(a, b) } -> std::same_as<T>;
    { dot_product(a, b) } -> std::same_as<scalar_type<T>>;
    { zero<T>() } -> std::same_as<T>;
};

template<typename T>
struct is_point : std::false_type {};

template<Point T>
struct is_point<T> : std::true_type {};

// Base case for grid type trait
template<Point T, std::size_t dim>
struct point_grid_type_trait {
    using type = typename point_grid_type_trait<std::vector<T>, dim - 1>::type;
};

// Specialization for dim = 1 (1D grid)
template<Point T>
struct point_grid_type_trait<T, 1> {
    using type = std::vector<T>;
};

// Specialization for Eigen::Matrix types (3D grid)
template<typename T, int Rows, int Cols, int Options, std::size_t dim>
struct point_grid_type_trait<Eigen::Matrix<T, Rows, Cols, Options>, dim> {
    using type = Eigen::Tensor<Eigen::Matrix<T, Rows, Cols, Options>, dim>;
};

// Helper alias to make syntax more straightforward when using the trait
template<typename T, std::size_t dim>
using point_grid_type = typename point_grid_type_trait<T, dim>::type;


// ================================================
// ============= Container Traits =================
// ================================================

// -------- Default Traits for std::vector --------

// Size Trait
template<typename T, typename = void>
struct size_trait {
    static std::size_t size(const T &container) { return container.size(); }
};

// Get Trait
template<typename T, typename = void>
struct get_trait {
    static auto get(const T &container, std::size_t i)
            -> decltype(container[i]) {
        return container[i];
    }
};

// Push Trait
template<typename T, typename = void>
struct push_trait {
    static void push(T &container, const typename T::value_type &value) {
        container.push_back(value);
    }
};

// Resize Trait
template<typename T, typename = void>
struct resize_trait {
    static void resize(T &container, std::size_t size) {
        container.resize(size);
    }
};

// ------- Specialization Traits for Eigen::Matrix -------

// Size Trait
template<typename T>
struct size_trait<
        T, std::enable_if_t<std::is_base_of_v<Eigen::MatrixBase<T>, T>>> {
    static std::size_t size(const T &container) { return container.cols(); }
};

// Get Trait
template<typename T>
struct get_trait<T,
                 std::enable_if_t<std::is_base_of_v<Eigen::MatrixBase<T>, T>>> {
    static auto get(const T &container, std::size_t i)
            -> decltype(container.col(i)) {
        return container.col(i);
    }
};

// Push Trait
template<typename U>
struct push_trait<Eigen::Matrix<U, 3, -1>> {
    static void push(Eigen::Matrix<U, 3, -1>      &container,
                     const Eigen::Matrix<U, 3, 1> &value, std::size_t i) {
        container.col(i) = value;
    }
};

// Resize Trait
template<typename T>
struct resize_trait<
        T, std::enable_if_t<std::is_base_of_v<Eigen::MatrixBase<T>, T>>> {
    static void resize(T &container, std::size_t size) {
        container.conservativeResize(Eigen::NoChange, size);
    }
};


//--- Eigen type specialisations ---//


//// Eigen squared_distance
//template<Point T>
//    requires std::derived_from<T, Eigen::EigenBase<T>>
//constexpr scalar_type<T> squared_distance(const T &a, const T &b) {
//    return (a - b).squaredNorm();
//}
//
//// Eigen squared_norm
//template<Point T>
//    requires std::derived_from<T, Eigen::EigenBase<T>>
//constexpr scalar_type<T> squared_norm(const T &a) {
//    return a.squaredNorm();
//}
//
//// Eigen element_wise_division
//template<Point T>
//    requires std::derived_from<T, Eigen::EigenBase<T>>
//constexpr T element_wise_division(const T &a, const T &b) {
//    return (a.array() / b.array()).matrix();
//}
//
//// Eigen zero
//template<Point T>
//    requires std::derived_from<T, Eigen::EigenBase<T>>
//constexpr T zero() {
//    return T::Zero();
//}
//
//
//// Eigen dot_product
//template<Point T>
//    requires std::derived_from<T, Eigen::EigenBase<T>>
//constexpr scalar_type<T> dot_product(const T &a, const T &b) {
//    return a.dot(b);
//}
//
//// Eigen normalize
//template<typename T>
//    requires std::same_as<T, Eigen::Matrix<scalar_type<T>, 3, 1>>
//constexpr Eigen::Matrix<scalar_type<T>, 3, 1> normalize(const T &a) {
//    return a.normalized();
//}
////--- Array type specialisations ---//
//
////    // Specialization for std::array<T, N>
////    template<typename T, std::size_t N>
////    struct point_series_type_trait<std::array<T, N>> {
////        using type = std::array<std::array<T, N>, N>; // Two-dimensional array
////    };
//
//
//// Array squared_distance
//template<Point T, std::size_t N>
//    requires std::same_as<T, std::array<scalar_type<T>, N>>
//constexpr scalar_type<T> squared_distance(const T &a, const T &b) {
//#ifdef __cpp_lib_transform_reduce
//    return std::transform_reduce(
//            a.begin(), a.end(), b.begin(), scalar_type<T>{0}, std::plus{},
//            [](auto a, auto b) { return std::pow(a - b, 2); });
//#else
//    scalar_type<T> sum = 0;
//    for (std::size_t i = 0; i < N; ++i) { sum += std::pow(a[i] - b[i], 2); }
//    return sum;
//#endif
//}
//
//// Array squared_norm
//template<Point T, std::size_t N>
//    requires std::same_as<T, std::array<scalar_type<T>, N>>
//constexpr scalar_type<T> squared_norm(const T &a) {
//#ifdef __cpp_lib_transform_reduce
//    return std::transform_reduce(a.begin(), a.end(), scalar_type<T>{0},
//                                 std::plus{},
//                                 [](auto a) { return std::pow(a, 2); });
//#else
//    scalar_type<T> sum = 0;
//    for (std::size_t i = 0; i < N; ++i) { sum += std::pow(a[i], 2); }
//    return sum;
//#endif
//}
//
//// Array element_wise_division
//template<Point T, std::size_t N>
//    requires std::same_as<T, std::array<scalar_type<T>, N>>
//constexpr T element_wise_division(const T &a, const T &b) {
//    T result;
//    std::ranges::transform(a, b, result.begin(), std::divides{});
//    return result;
//}
//
//template<Point T, std::size_t N>
//    requires std::same_as<T, std::array<scalar_type<T>, N>>
//constexpr scalar_type<T> dot_product(const T &a, const T &b) {
//#ifdef __cpp_lib_transform_reduce
//    return std::transform_reduce(a.begin(), a.end(), b.begin(),
//                                 scalar_type<T>{0}, std::plus{},
//                                 std::multiplies{});
//#else
//    scalar_type<T> sum = 0;
//    for (std::size_t i = 0; i < N; ++i) { sum += a[i] * b[i]; }
//    return sum;
//#endif
//}
//
//// Array zero
//template<Point T, std::size_t N>
//    requires std::same_as<T, std::array<scalar_type<T>, N>>
//constexpr T zero() {
//    return T{};
//}
//
//// Array Euclidean norm
//template<typename T, std::size_t N>
//    requires std::same_as<T, std::array<scalar_type<T>, N>>
//constexpr scalar_type<T> normalize(const T &a) {
//    return std::sqrt(std::accumulate(a.begin(), a.end(),
//                                     static_cast<scalar_type<T>>(0),
//                                     [](const auto &sum, const auto &element) {
//                                         return sum + element * element;
//                                     }));
//}

}// namespace odin::shape

#endif//ODIN_SHAPE_CONCEPTS_H
