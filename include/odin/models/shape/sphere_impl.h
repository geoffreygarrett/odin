#ifndef ODIN_SPHERE_IMPL_H
#define ODIN_SPHERE_IMPL_H

#include "shape_concepts.h"
#include "sphere.h"
#include <cmath>
#include <odin/core/point.h>

namespace odin::shape {

// constructor
template<typename U, Point P>
sphere<U, P>::sphere(U radius, P center) : base(center), m_radius(radius) {}

template<typename U, Point P>
U sphere<U, P>::volume_impl() {
    return (4.0 / 3.0) * M_PI * std::pow(m_radius, 3);
}

template<typename U, Point P>
U sphere<U, P>::surface_area_impl() {
    return 4 * M_PI * std::pow(m_radius, 2);
}

template<typename U, Point P>
P sphere<U, P>::centroid_impl() {
    return m_center;
}

template<typename U, Point P>
bool sphere<U, P>::is_inside_impl(const P &point) {
    U distance_squared = squared_distance(point, m_center);
    return distance_squared <= std::pow(m_radius, 2);
}

// sphere specific functions
template<typename U, Point P>
U sphere<U, P>::get_radius() const {
    return m_radius;
}

template<typename U, Point P>
void sphere<U, P>::set_radius(U radius) {
    m_radius = radius;
}

template<typename U, Point P>
//    std::optional<P> sphere<U, P>::ray_intersection_impl(const Ray<P> &ray) {
std::optional<P> sphere<U, P>::ray_intersection_impl(const P &origin,
                                                     const P &direction) {
    // Vector from the sphere center to the ray origin
    P oc = origin - m_center;
    U a  = squared_norm(direction);
    U b  = 2.0 * dot_product(oc, direction);
    U c  = squared_norm(oc) - m_radius * m_radius;

    U discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return std::nullopt;// No intersection
    } else {
        U t = (-b - std::sqrt(discriminant)) / (2.0 * a);
        if (t < 0) {
            t = (-b + std::sqrt(discriminant))
              / (2.0
                 * a);// If the closer intersection point is behind the ray, choose the farther one
        }
        return origin + t * direction;// Intersection point
    }
}
//    template<typename U, Point P>
//    // Returns the intersection point, or NaN if there is no intersection.
//    __device__ P sphere<U, P>::ray_intersection_impl(const P &origin, const P &direction) {
//        // Vector from the sphere center to the ray origin
//        P oc = origin - m_center;
//        U a  = squared_norm(direction);
//        U b  = 2.0 * dot_product(oc, direction);
//        U c  = squared_norm(oc) - m_radius * m_radius;
//
//        U discriminant = b * b - 4 * a * c;
//        if (discriminant < 0) {
//            return P(std::numeric_limits<U>::quiet_NaN()); // No intersection
//        } else {
//            U t = (-b - std::sqrt(discriminant)) / (2.0 * a);
//            if (t < 0) {
//                t = (-b + std::sqrt(discriminant)) / (2.0 * a);
//                // If the closer intersection point is behind the ray, choose the farther one
//            }
//            return origin + t * direction; // Intersection point
//        }
//    }

}// namespace odin::shape

#endif//ODIN_SPHERE_IMPL_H
