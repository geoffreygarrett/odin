#ifndef ODIN_FACE_VERTEX_H
#define ODIN_FACE_VERTEX_H

#include "mesh_base.h"
#include <vector>

namespace odin::mesh {

template<typename U, typename P, std::size_t dim = 3>
class face_vertex : public mesh_base<face_vertex<U, P>, U, P, dim> {
public:
    using base        = mesh_base<face_vertex<U, P>, U, P, dim>;
    using scalar_type = typename base::scalar_type;
    using point_type  = typename base::point_type;
    using facet_type  = typename base::facet_type;

    std::vector<point_type> vertices;
    std::vector<facet_type> faces;

    explicit face_vertex(point_type center = zero<point_type>());

    facet_type get_facet(std::size_t index) const;

    std::vector<point_type> get_vertices() const;

    std::vector<facet_type> get_facets() const;

    [[nodiscard]] std::size_t get_face_count() const;
};
}// namespace odin::mesh

#include "face_vertex_impl.h"

#endif//ODIN_FACE_VERTEX_H
