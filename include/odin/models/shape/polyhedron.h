// polyhedron.h

#ifndef ODIN_POLYHEDRON_H
#define ODIN_POLYHEDRON_H

#include "shape_base.h"
#include "shape_concepts.h"

namespace odin::shape {
    template<typename U, Point P>
    class polyhedron : public shape_base<polyhedron<U, P>, U, P, 3> {
    public:
        // base
        using base = shape_base<polyhedron<U, P>, U, P, 3>;

        using point_type  = P;
        using index_array = std::vector<std::size_t>;
        using facet_type  = std::vector<point_type>;

        // constructor
        explicit polyhedron(U radius_a, U radius_b, U radius_c,
                            P center = zero<P>());
    };
}// namespace odin::shape

#endif// ODIN_POLYHEDRON_H