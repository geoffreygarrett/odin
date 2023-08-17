#ifndef ODIN_POLYHEDRAL_HPP
#define ODIN_POLYHEDRAL_HPP

#include "gravitational_base.hpp"
#include <polyhedralGravity/calculation/GravityModel.h>
#include <polyhedralGravity/calculation/MeshChecking.h>
#include <polyhedralGravity/input/TetgenAdapter.h>
#include <polyhedralGravity/model/GravityModelData.h>
#include <polyhedralGravity/model/Polyhedron.h>


namespace odin {
    using namespace polyhedralGravity;

    template<typename Scalar, std::size_t Dim = 3>
    class Polyhedral : public GravitationalModelBase<Polyhedral<Scalar, Dim>, Scalar, Dim> {
    public:
        using scalar_type     = Scalar;
        using vector_type     = Eigen::Matrix<Scalar, Dim, 1>;
        using matrix_type     = Eigen::Matrix<Scalar, Dim, Dim>;
        using polyhedron_type = std::unique_ptr<Polyhedron>;
        using nodes_type      = std::vector<std::array<double, 3>>;
        using faces_type      = std::vector<std::array<size_t, 3>>;

        Polyhedral() = default;

        Polyhedral(const Polyhedral &other)
            : m_density(other.m_density),
              m_polyhedron(std::make_unique<Polyhedron>(*other.m_polyhedron)) {}

        //        Polyhedral(std::unique_ptr<Polyhedron> polyhedron, const Scalar &density)
        //            : m_polyhedron(std::move(polyhedron)), m_density(density) {}

        Polyhedral thread_local_copy_impl() const {
            return Polyhedral(*this);
        }

        Polyhedral(const std::string &file_path, const Scalar &density)
            : m_density(density) {
            TetgenAdapter adapter({file_path});
            m_polyhedron = std::make_unique<Polyhedron>(adapter.getPolyhedron());
        }

        Polyhedral(nodes_type nodes, faces_type faces, scalar_type density)
            : m_density(density) {
            m_polyhedron = std::make_unique<Polyhedron>(nodes, faces);
        }

        Polyhedral(const std::vector<std::string> &file_paths, const Scalar &density)
            : m_density(density) {
            TetgenAdapter adapter(file_paths);
            m_polyhedron = std::make_unique<Polyhedron>(adapter.getPolyhedron());
        }

        scalar_type potential_impl(const vector_type &position) const {
            auto result = evaluate_gravity_model(position);
            return static_cast<Scalar>(result.gravitationalPotential);
        }

        vector_type acceleration_impl(const vector_type &position) const {
            auto result = evaluate_gravity_model(position);
            // TODO: Plots were not correct, perhaps a sign error?
            //  Check if this is my mistake, or the library's mistake.
            return Eigen::Map<const vector_type>(result.acceleration.data());
        }

        matrix_type tensor_impl(const vector_type &position) {
            auto result = evaluate_gravity_model(position);
            return convert_to_matrix(result.gradiometricTensor);
        }

        std::tuple<scalar_type, vector_type, matrix_type> evaluate_imp(const vector_type &position) const {
            auto result = evaluate_gravity_model(position);
            return std::make_tuple(
                    static_cast<Scalar>(result.gravitationalPotential),
                    Eigen::Map<const vector_type>(result.acceleration.data()),
                    convert_to_matrix(result.gradiometricTensor));
        }

    private:
        polyhedron_type m_polyhedron;
        scalar_type     m_density;

        GravityModelResult evaluate_gravity_model(const vector_type &position) const {
            std::array<double, 3> computation_point{};// Convert from vector_type to Array3
            computation_point[0] = position[0];
            computation_point[1] = position[1];
            computation_point[2] = position[2];
            return polyhedralGravity::GravityModel::evaluate(*m_polyhedron, static_cast<double>(m_density), computation_point);
        }

        matrix_type convert_to_matrix(const std::array<double, 6> &gradiometric_tensor) {
            matrix_type tensor;
            // Convert std::array<double, 6> to matrix_type. Depending on the layout of the data, this might need adjustment.
            tensor(0, 0) = gradiometric_tensor[0];
            tensor(0, 1) = gradiometric_tensor[1];
            tensor(0, 2) = gradiometric_tensor[2];
            tensor(1, 0) = gradiometric_tensor[3];
            tensor(1, 1) = gradiometric_tensor[4];
            tensor(1, 2) = gradiometric_tensor[5];
            return tensor;
        }
    };
}// namespace odin

#endif// ODIN_POLYHEDRAL_HPP