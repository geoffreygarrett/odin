#ifndef GRAVITATIONAL_BASE_HPP
#define GRAVITATIONAL_BASE_HPP

#include <tbb/tbb.h>

#include <Eigen/Core>
#include <odin/parallel/parallel.hpp>
#include <odin/parallel/parallel_concepts.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

#include "gravitational_concept.hpp"

namespace odin {

//#include <autodiff/forward.hpp>

template<typename Derived>
class crtp_base {
public:
    Derived &as_derived() {
        return static_cast<Derived &>(*this);
    }

    const Derived &as_derived() const {
        return static_cast<const Derived &>(*this);
    }
};

// autodiff policy
//template<typename Derived, typename Scalar, std::size_t Dim = 3>

template<typename Derived, typename Scalar, std::size_t Dim = 3>
class GravitationalModelBase : public crtp_base<Derived> {
public:
    using vector_type = Eigen::Matrix<Scalar, Dim, 1>;
    using crtp_base<Derived>::as_derived;

    Derived thread_local_copy() {
        return as_derived().thread_local_copy_impl();
    }

    Scalar potential(const vector_type &position) const {
        return as_derived().potential_impl(position);
    }

    // enalbe if 'Scalar' is type autodiff::dual
    //    template<typename T = Scalar>
    //        requires std::is_same_v<T, autodiff::dual>
    //    vector_type partial_potential_wrt_position(
    //            const vector_type &position, const vector_type &position2) const {
    //        return as_derived().partial_potential_wrt_position_impl(position, position2);
    //    }

    vector_type acceleration(const vector_type &position) const {
        return as_derived().acceleration_impl(position);
    }

    // enalbe if 'Scalar' is type autodiff::dual
    //    template<typename T = Scalar>
    //        requires std::is_same_v<T, autodiff::dual>
    //    vector_type partial_acceleration_wrt_position(
    //            const vector_type &position, const vector_type &position2) const {
    //        return as_derived().partial_acceleration_wrt_position_impl(position, position2);
    //    }

    using series_vector_type = Eigen::Matrix<Scalar, Eigen::Dynamic, Dim>;
    using series_scalar_type = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // using series_vector_type = Eigen::Matrix<Scalar, Eigen::Dynamic, Dim>;
    // using series_scalar_type = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // Series version, potential
    series_scalar_type potential_series(const series_vector_type &positions)
        requires ThreadSafe<Derived>
    {
        spdlog::debug("Calculating potential in series");
        series_scalar_type potentials(positions.rows());

        SeriesCalculator<vector_type,// serial input
                Scalar,              // serial output
                series_vector_type,  // parallel input
                series_scalar_type   // parallel output
                >::
                calculate(
                        as_derived(),
                        positions,
                        [](Derived &model, const vector_type &v) -> Scalar {
                            return model.potential(v);
                        },
                        potentials);

        return potentials;
    }

    // Series version, acceleration
    series_vector_type acceleration_series(const series_vector_type &positions)
        requires ThreadSafe<Derived>
    {
        spdlog::debug("Calculating acceleration in series");
        series_vector_type accelerations(positions.rows(), positions.cols());

        // Call SeriesCalculator with input_dims and output_dims
        SeriesCalculator<vector_type,// serial input
                vector_type,         // serial output
                series_vector_type,  // parallel input
                series_vector_type   // parallel output
                >::
                calculate(
                        as_derived(),
                        positions,
                        [](Derived &model, const vector_type &v) -> vector_type {
                            return model.acceleration(v);
                        },
                        accelerations);

        return accelerations;
    }

    using grid_3d_scalar_type = Eigen::Tensor<Scalar, 3>;
    grid_3d_scalar_type potential_grid(const std::array<grid_3d_scalar_type, 3> &grid_positions)
        requires ThreadSafe<Derived>
    {
        grid_3d_scalar_type potential(grid_positions[0].dimension(0),
                grid_positions[1].dimension(1),
                grid_positions[2].dimension(2));

        GridCalculator<3, 1>::calculate(
                as_derived(),
                grid_positions,
                [](Derived &model, const vector_type &v) -> Scalar { return model.potential(v); },
                potential);

        return potential;
    }

    using grid_2d_scalar_type = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    grid_2d_scalar_type potential_grid(const std::array<grid_2d_scalar_type, 3> &grid_positions)
        requires ThreadSafe<Derived>
    {
        grid_2d_scalar_type potential(grid_positions[0].rows(), grid_positions[0].cols());

        GridCalculator<2, 1>::calculate(
                as_derived(),
                grid_positions,
                [](Derived &model, const vector_type &v) -> Scalar { return model.potential(v); },
                potential);

        return potential;
    }

    std::array<grid_3d_scalar_type, 3> acceleration_grid(
            const std::array<grid_3d_scalar_type, 3> &grid_positions)
        requires ThreadSafe<Derived>
    {
        std::array<grid_3d_scalar_type, 3> accelerations;
        for (auto &a: accelerations) {
            a = grid_3d_scalar_type(grid_positions[0].dimension(0),
                    grid_positions[1].dimension(1),
                    grid_positions[2].dimension(2));
        }

        GridCalculator<3, 3>::calculate(
                as_derived(),
                grid_positions,
                [](Derived &model, vector_type &v) -> vector_type { return model.acceleration(v); },
                accelerations);

        return accelerations;
    }
};

}// namespace odin
#endif// GRAVITATIONAL_BASE_HPP
