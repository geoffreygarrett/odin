#ifndef ODIN_MESH_BASE_H
#define ODIN_MESH_BASE_H

#include <array>
#include <odin/core/point.h>
#include <vector>

namespace odin::mesh {
using namespace odin::shape;

template<typename T, typename U, typename P, std::size_t dim>
class mesh_base {
public:
    using scalar_type   = U;
    using point_type    = P;// User-provided point container
    using index_type    = std::size_t;
    using vertex_type   = std::array<point_type, dim>;
    using facet_type    = std::array<index_type, dim>;
    using vertices_type = std::vector<vertex_type>;
    using facets_type   = std::vector<facet_type>;

    // constructor
    explicit mesh_base(point_type center = zero<point_type>())
        : m_center(std::move(center)) {}

private:
    point_type m_center;
};

}// namespace odin::mesh

#endif//ODIN_MESH_BASE_H
