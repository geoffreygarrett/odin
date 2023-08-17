// ellipsoid_impl.h

#ifndef ODIN_ELLIPSOID_IMPL_H
#define ODIN_ELLIPSOID_IMPL_H

#include "ellipsoid.h"
#include "shape_concepts.h"
#include <cmath>

namespace odin::shape {

template<typename U, Point P>
ellipsoid<U, P>::ellipsoid(U radius_a, U radius_b, U radius_c, P center)
    : m_radius_a(radius_a), m_radius_b(radius_b), m_radius_c(radius_c),
      base(center) {}

template<typename U, Point P>
ellipsoid<U, P>::ellipsoid(P radii, P center)
    : m_radius_a(radii[0]), m_radius_b(radii[1]), m_radius_c(radii[2]),
      base(center) {}

template<typename U, Point P>
U ellipsoid<U, P>::volume_impl() {
    return (4.0 / 3.0) * M_PI * m_radius_a * m_radius_b * m_radius_c;
}

template<typename U, Point P>
U ellipsoid<U, P>::surface_area_impl() {
    constexpr double p = 1.6075;
    return 4 * M_PI
         * std::cbrt((std::pow(m_radius_a, p) * std::pow(m_radius_b, p)
                      + std::pow(m_radius_a, p) * std::pow(m_radius_c, p)
                      + std::pow(m_radius_b, p) * std::pow(m_radius_c, p))
                     / 3.0);
}

template<typename U, Point P>
bool ellipsoid<U, P>::is_inside_impl(const P &point) {
    P relative_point   = point - m_center;
    P normalized_point = element_wise_division(
            relative_point, P(m_radius_a, m_radius_b, m_radius_c));
    U distance_normalized = squared_norm(normalized_point);
    return distance_normalized <= 1;
}

// Get silhouette dimensions (semi-major and semi-minor axes)
// from the perspective of a viewer located at 'source' looking towards 'direction'.
template<typename U, Point P>
std::pair<U, U>
ellipsoid<U, P>::get_silhouette_dimensions(const P &source,
                                           const P &direction) const {
    // Shift the world space to the ellipsoid's local space
    P localSource    = source - m_center;
    P localDirection = direction - m_center;

    // Compute the view direction
    P view_direction_vec = localDirection - localSource;
    P view_direction     = normalize(view_direction_vec);

    // Compute the length of the ellipsoid's radii in the view direction
    U length_a
            = dot_product(P(m_radius_a * P{1, 0, 0}), view_direction);
    U length_b
            = dot_product(P(m_radius_b * P{0, 1, 0}), view_direction);
    U length_c
            = dot_product(P(m_radius_c * P{0, 0, 1}), view_direction);

    // The lengths form a right triangle with the direction vector, so we can use the Pythagorean theorem to find the silhouette dimensions.
    // The semi-major axis (a_prime) will be the hypotenuse of this triangle, and the semi-minor axis (b_prime) will be the other side.
    U a_prime = std::sqrt(length_a * length_a + length_b * length_b);
    U b_prime = length_c;

    return {a_prime, b_prime};
}

template<typename U, Point P>
P ellipsoid<U, P>::centroid_impl() {
    return m_center;
}

// Ellipsoid specific functions
template<typename U, Point P>
U ellipsoid<U, P>::get_radius_a() const {
    return m_radius_a;
}

template<typename U, Point P>
void ellipsoid<U, P>::set_radius_a(U radius_a) {
    m_radius_a = radius_a;
}

template<typename U, Point P>
U ellipsoid<U, P>::get_radius_b() const {
    return m_radius_b;
}

template<typename U, Point P>
void ellipsoid<U, P>::set_radius_b(U radius_b) {
    m_radius_b = radius_b;
}

template<typename U, Point P>
U ellipsoid<U, P>::get_radius_c() const {
    return m_radius_c;
}

template<typename U, Point P>
void ellipsoid<U, P>::set_radius_c(U radius_c) {
    m_radius_c = radius_c;
}

template<typename U, Point P>
P ellipsoid<U, P>::get_radii() const {
    return P(m_radius_a, m_radius_b, m_radius_c);
}

template<typename U, Point P>
void ellipsoid<U, P>::set_radii(P radii) {
    m_radius_a = radii[0];
    m_radius_b = radii[1];
    m_radius_c = radii[2];
}

template<typename U, Point P>
//    std::optional<P> ellipsoid<U, P>::ray_intersect_impl(const Ray<P> &ray) {
std::optional<P>
ellipsoid<U, P>::ray_intersection_impl(const P &ray_origin,
                                       const P &ray_direction) {
    // Direction vector adjusted for the ellipsoid's radii
    P dir = element_wise_division(ray_direction, get_radii());

    // Vector from the ellipsoid center to the ray origin, adjusted for the ellipsoid's radii
    P oc = element_wise_division(P(ray_origin - m_center), get_radii());

    U a = squared_norm(dir);
    U b = 2.0 * dot_product(oc, dir);
    U c = squared_norm(oc) - 1;

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
        return ray_origin + t * ray_direction;// Intersection point
    }
}


//template<typename U, typename P>
//#if defined(__CUDACC__)
//__host__ __device__
//#endif
//        thrust::pair<P, bool>
//        ellipsoid<U, P>::ray_intersection_impl(const P& ray_origin,
//                                               const P& ray_direction,
//                                               const P& ellipsoid_center,
//                                               const P& ellipsoid_radii)
//{
//    // Direction vector adjusted for the ellipsoid's radii
//    P dir = element_wise_division(ray_direction, ellipsoid_radii);
//
//    // Vector from the ellipsoid center to the ray origin, adjusted for the ellipsoid's radii
//    P oc = element_wise_division(ray_origin - ellipsoid_center, ellipsoid_radii);
//
//    U a = squared_norm(dir);
//    U b = 2.0 * dot_product(oc, dir);
//    U c = squared_norm(oc) - 1;
//
//    U discriminant = b * b - 4 * a * c;
//    if (discriminant < 0) {
//        return thrust::make_pair(P(), false); // No intersection
//    }
//    else {
//        U t = (-b - std::sqrt(discriminant)) / (2.0 * a);
//        if (t < 0) {
//            t = (-b + std::sqrt(discriminant)) / (2.0 * a); // If the closer intersection point is behind the ray, choose the farther one
//        }
//        return thrust::make_pair(ray_origin + t * ray_direction, true); // Intersection point
//    }
//}


}// namespace odin::shape

#endif//ODIN_ELLIPSOID_IMPL_H
