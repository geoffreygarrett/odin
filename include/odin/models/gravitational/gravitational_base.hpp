#ifndef GRAVITATIONAL_BASE_HPP
#define GRAVITATIONAL_BASE_HPP

#include <Eigen/Core>

template<typename T, size_t Dim = 3>
concept GravitationalConcept = requires(T a, const Eigen::Vector<typename T::Scalar, Dim> &position) {
    // Ensure that a type T has a type member Scalar
    typename T::Scalar;

    // Ensure that a member function potential() exists with the appropriate signature
    { a.potential(position) } -> std::convertible_to<typename T::Scalar>;

    // Ensure that a member function acceleration() exists with the appropriate signature
    { a.acceleration(position) } -> std::convertible_to<Eigen::Vector<typename T::Scalar, Dim>>;
};

// Helper function to ensure gravitational models are compatible with the GravitationalModel concept
template<typename T, size_t Dim = 3>
constexpr bool is_gravitational_model_v = GravitationalConcept<T, Dim>;

#endif// GRAVITATIONAL_BASE_HPP
