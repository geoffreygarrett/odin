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

//--- Array type specialisations ---//

}// namespace odin::traits

#endif//ODIN_SHAPE_CONCEPTS_H
