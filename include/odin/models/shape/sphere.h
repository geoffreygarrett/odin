#ifndef ODIN_SPHERE_H
#define ODIN_SPHERE_H

#include "shape_base.h"
#include "shape_concepts.h"// include the utility namespace
#include <type_traits>

namespace odin::shape {

    template<typename U, Point P>
    class sphere : public shape_base<sphere<U, P>, U, P, 3> {
    public:
        // base class
        using base = shape_base<sphere<U, P>, U, P, 3>;

        // constructor
        explicit sphere(U radius, P center = zero<P>());

        // interface functions
        U    volume_impl();
        U    surface_area_impl();
        P    centroid_impl();
//        std::optional<P> ray_intersection_impl(const Ray<P> &ray);
        std::optional<P> ray_intersection_impl(const P &origin, const P &direction);
        bool is_inside_impl(const P &point);

        // sphere specific functions
        U    get_radius() const;
        void set_radius(U radius);


    private:
        U m_radius;
        P m_center;
    };

}// namespace odin


#endif//ODIN_SPHERE_H
