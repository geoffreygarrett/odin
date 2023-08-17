// ellipsoid.h

#ifndef ODIN_ELLIPSOID_H
#define ODIN_ELLIPSOID_H

#include "shape_base.h"
#include "shape_concepts.h"

namespace odin::shape {

template<typename U, Point P>
class ellipsoid : public shape_base<ellipsoid<U, P>, U, P, 3> {
public:
    // base
    using base = shape_base<ellipsoid<U, P>, U, P, 3>;

    // Constructor
    explicit ellipsoid(
            U radius_a, U radius_b, U radius_c, P center = zero<P>());

    explicit ellipsoid(P radii, P center = zero<P>());

    // Implementation of the interface
    U                volume_impl();
    U                surface_area_impl();
    P                centroid_impl();
    bool             is_inside_impl(const P &point);
    std::optional<P> ray_intersection_impl(const P &origin, const P &direction);


    std::pair<U, U> get_silhouette_dimensions(
            const P &source, const P &direction) const;

    // Ellipsoid specific functions
    U    get_radius_a() const;
    void set_radius_a(U radius_a);
    U    get_radius_b() const;
    void set_radius_b(U radius_b);
    U    get_radius_c() const;
    void set_radius_c(U radius_c);
    P    get_radii() const;
    void set_radii(P radii);


private:
    U m_radius_a;
    U m_radius_b;
    U m_radius_c;
    P m_center;
};

}// namespace odin::shape

#endif//ODIN_ELLIPSOID_H