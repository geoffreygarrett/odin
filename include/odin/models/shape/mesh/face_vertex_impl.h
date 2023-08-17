#ifndef ODIN_FACE_VERTEX_IMPL_H
#define ODIN_FACE_VERTEX_IMPL_H

#include "face_vertex.h"

namespace odin::mesh {

template<typename U, typename P, std::size_t dim>
face_vertex<U, P, dim>::face_vertex(point_type center) : base(center) {}

template<typename U, typename P, std::size_t dim>
typename face_vertex<U, P, dim>::facet_type face_vertex<U, P, dim>::get_facet(
        std::size_t index) const {
    return faces[index];
}

template<typename U, typename P, std::size_t dim>
std::size_t face_vertex<U, P, dim>::get_face_count() const {
    return faces.size();
}

template<typename U, typename P, std::size_t dim>
std::vector<typename face_vertex<U, P, dim>::point_type>
face_vertex<U, P, dim>::get_vertices() const {
    return vertices;
}

template<typename U, typename P, std::size_t dim>
std::vector<typename face_vertex<U, P, dim>::facet_type>
face_vertex<U, P, dim>::get_facets() const {
    return faces;

}
}// namespace odin::mesh

#endif//ODIN_FACE_VERTEX_IMPL_H