#ifndef GRAVITATIONAL_BASE_HPP
#define GRAVITATIONAL_BASE_HPP

#include <Eigen/Core>
#include <tbb/tbb.h>
#include <unsupported/Eigen/CXX11/Tensor>


template<typename T, size_t Dim = 3>
concept GravitationalConcept = requires(T a, const Eigen::Vector<typename T::Scalar, Dim> &position) {
    // Ensure that a type T has a type member Scalar
    typename T::Scalar;

    // Ensure that a member function potential() exists with the appropriate signature
    { a.potential(position) } -> std::convertible_to<typename T::Scalar>;

    // Ensure that a member function acceleration() exists with the appropriate signature
    { a.acceleration(position) } -> std::convertible_to<Eigen::Vector<typename T::Scalar, Dim>>;
};

// Helper function to ensure gravitational models are compatible with the GravitationalModel concept
template<typename T, size_t Dim = 3>
constexpr bool is_gravitational_model_v = GravitationalConcept<T, Dim>;


template<typename Derived, typename Scalar, int Dim = 3>
class GravitationalModelBase {
public:
    using VectorType = Eigen::Matrix<Scalar, Dim, 1>;

// destructor
    virtual ~GravitationalModelBase() = default;

    // 3D version
    Eigen::Tensor<Scalar, 3> calculate_potentials(
            const Eigen::Tensor<Scalar, 3> &x_grid,
            const Eigen::Tensor<Scalar, 3> &y_grid,
            const Eigen::Tensor<Scalar, 3> &z_grid
    ) const {
        assert(x_grid.dimensions() == y_grid.dimensions() && x_grid.dimensions() == z_grid.dimensions());

        const auto &dims = x_grid.dimensions();

        Eigen::Tensor<Scalar, 3> potential_grid(dims);

        // Use lambda function to wrap member function
        auto fn = [this](Derived &model, const VectorType &v) -> Scalar {
            return (model.*(&Derived::potential))(v);
        };

        // parallel calculation
        calculate_in_parallel<Eigen::Tensor<Scalar, 3>>(
                x_grid, y_grid, z_grid, fn, potential_grid
        );

        return potential_grid;
    }

    // 3D version
    Eigen::Tensor<Scalar, 4> calculate_accelerations(
            const Eigen::Tensor<Scalar, 3> &x_grid,
            const Eigen::Tensor<Scalar, 3> &y_grid,
            const Eigen::Tensor<Scalar, 3> &z_grid
    ) const {
        assert(x_grid.dimensions() == y_grid.dimensions() && x_grid.dimensions() == z_grid.dimensions());

        const auto &dims = x_grid.dimensions();

        Eigen::Tensor<Scalar, 4> acceleration_grid(dims[0], dims[1], dims[2], Dim);

        // Use lambda function to wrap member function
        auto fn = [this](Derived &model, const VectorType &v) -> Eigen::Matrix<Scalar, Dim, 1> {
            return (model.*(&Derived::acceleration))(v);
        };

        // parallel calculation
        calculate_in_parallel<Eigen::Tensor<Scalar, 4>>(
                x_grid, y_grid, z_grid, fn, acceleration_grid
        );

        return acceleration_grid;
    }

    Derived thread_local_copy() const {
        return static_cast<const Derived &>(*this).thread_local_copy_impl();
    }


private:
    template<typename OutputTensorType, typename Callable>
    void calculate_in_parallel(
            const Eigen::Tensor<Scalar, 3> &x_grid,
            const Eigen::Tensor<Scalar, 3> &y_grid,
            const Eigen::Tensor<Scalar, 3> &z_grid,
            Callable &&calculate_fn,
            OutputTensorType &output_grid
    ) const {
        static_assert(OutputTensorType::NumDimensions == 3 || OutputTensorType::NumDimensions == 4,
                      "OutputTensorType must be 3D or 4D.");

        const auto &dims = x_grid.dimensions();

        tbb::parallel_for(
                tbb::blocked_range3d<int>(0, dims[0], 0, dims[1], 0, dims[2]),
                [&](const tbb::blocked_range3d<int> &r) {
                    thread_local Derived tl_model = static_cast<const Derived &>(*this).thread_local_copy();

                    for (int i = r.pages().begin(); i != r.pages().end(); ++i) {
                        for (int j = r.rows().begin(); j != r.rows().end(); ++j) {
                            for (int k = r.cols().begin(); k != r.cols().end(); ++k) {
                                VectorType position(x_grid(i, j, k), y_grid(i, j, k), z_grid(i, j, k));

                                if constexpr (OutputTensorType::NumDimensions == 4) {
                                    // output tensor for acceleration
                                    Eigen::Matrix<Scalar, Dim, 1> result = calculate_fn(tl_model, position);
                                    for (int dim = 0; dim < Dim; ++dim) {
                                        output_grid(i, j, k, dim) = result[dim];
                                    }
                                } else if constexpr (OutputTensorType::NumDimensions == 3) {
                                    // output tensor for potential
                                    output_grid(i, j, k) = calculate_fn(tl_model, position);
                                }
                            }
                        }
                    }
                }
        );
    }

};

#endif// GRAVITATIONAL_BASE_HPP
