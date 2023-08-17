#ifndef ODIN_EIGEN_POINT_IMPL_H
#define ODIN_EIGEN_POINT_IMPL_H

#include <odin/core/point.h>
#include <type_traits>

namespace odin::shape {

// Eigen squared_distance
template<Point T>
    requires std::derived_from<T, Eigen::EigenBase<T>>
constexpr scalar_type<T> squared_distance(const T &a, const T &b) {
    return (a - b).squaredNorm();
}

// Eigen squared_norm
template<Point T>
    requires std::derived_from<T, Eigen::EigenBase<T>>
constexpr scalar_type<T> squared_norm(const T &a) {
    return a.squaredNorm();
}

// Eigen element_wise_division
template<Point T>
    requires std::derived_from<T, Eigen::EigenBase<T>>
constexpr T element_wise_division(const T &a, const T &b) {
    return (a.array() / b.array()).matrix();
}

// Eigen zero
template<Point T>
    requires std::derived_from<T, Eigen::EigenBase<T>>
constexpr T zero() {
    return T::Zero();
}

// Eigen dot_product
template<Point T>
    requires std::derived_from<T, Eigen::EigenBase<T>>
constexpr scalar_type<T> dot_product(const T &a, const T &b) {
    return a.dot(b);
}

// Eigen normalize
template<typename T>
    requires std::same_as<T, Eigen::Matrix<scalar_type<T>, 3, 1>>
constexpr Eigen::Matrix<scalar_type<T>, 3, 1> normalize(const T &a) {
    return a.normalized();
}

}// namespace odin::shape

#endif//ODIN_EIGEN_POINT_IMPL_H