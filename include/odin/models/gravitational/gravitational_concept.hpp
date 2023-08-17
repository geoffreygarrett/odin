#ifndef ODIN_GRAVIATIONAL_CONCEPT_H
#define ODIN_GRAVIATIONAL_CONCEPT_H

#include <Eigen/Core>
#include <cstdlib>

namespace odin {

    template<typename T, std::size_t Dim = 3>
    concept GravitationalConcept = requires(T a, const Eigen::Vector<typename T::Scalar, Dim> &position) {
        // Ensure that a type T has a type member Scalar
        typename T::Scalar;

        // Ensure that a member function potential() exists with the appropriate signature
        { a.potential(position) } -> std::convertible_to<typename T::Scalar>;

        // Ensure that a member function acceleration() exists with the appropriate signature
        { a.acceleration(position) } -> std::convertible_to<Eigen::Vector<typename T::Scalar, Dim>>;

        //    // Ensure that a member function spacial_gradient() exists with the appropriate signature
        //    { a.gradient(position) } -> std::convertible_to<Eigen::Matrix<typename T::Scalar, Dim, Dim>>;
    };

    // Helper function to ensure gravitational models are compatible with the GravitationalModel concept
    template<typename T, size_t Dim = 3>
    constexpr bool is_gravitational_model_v = GravitationalConcept<T, Dim>;

}// namespace odin


#endif//ODIN_GRAVIATIONAL_CONCEPT_H

//Simple: Point Mass​
//
//                 Simple: Hollow Sphere​
//
//                                  Simple: Homogeneous Sphere​
//
//                                                   Tri-axial Ellipsoid​
//
//                                                   Spherical Harmonic​
//
//                                                           Polyhedral