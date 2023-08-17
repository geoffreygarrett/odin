#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <odin/logging.hpp>
#include <odin/parallel/grid.hpp>
#include <odin/parallel/series.hpp>
#include <unsupported/Eigen/CXX11/Tensor>


namespace {
    using namespace odin;
    struct MockModel {
        int value;

        explicit MockModel(int v) : value(v) {}

        using scalar_type        = int;
        static constexpr int Dim = 3;
        using vector_type        = Eigen::Vector<scalar_type, Dim>;
        using series_vector_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, Dim>;
        using series_scalar_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

        // thread local copy
        [[nodiscard]] MockModel thread_local_copy() const {
            return MockModel(value);
        }

        [[nodiscard]] scalar_type potential(const vector_type &v) const {
            // Dummy implementation, replace with your logic
            return value * v.sum();
        }

        [[nodiscard]] vector_type acceleration(const vector_type &v) const {
            // Dummy implementation, replace with your logic
            return vector_type::Constant(value * v.sum());
        }

        // Series version
        series_scalar_type potential_series(
                const series_vector_type &positions) {
            spdlog::info("Calculating potential in series");

            series_scalar_type potentials(positions.rows());

            // Call SeriesCalculator with input_dims and 1 (for scalar output)
            SeriesCalculator<
                    vector_type,       // serial input
                    scalar_type,       // serial output
                    series_vector_type,// parallel input
                    series_scalar_type // parallel output
                    >::calculate(
                    *this, positions, [](MockModel &model, const vector_type &v) -> scalar_type {
                        return model.potential(v);
                    },
                    potentials);

            return potentials;
        }

        // Series version
        series_vector_type acceleration_series(
                const series_vector_type &positions) {
            spdlog::info("Calculating acceleration in series");
            series_vector_type accelerations(positions.rows(), positions.cols());

            // Call SeriesCalculator with input_dims and output_dims
            SeriesCalculator<
                    vector_type,       // serial input
                    vector_type,       // serial output
                    series_vector_type,// parallel input
                    series_vector_type // parallel output
                    >::calculate(
                    *this, positions, [](MockModel &model, const vector_type &v) -> vector_type {
                        return model.acceleration(v);
                    },
                    accelerations);

            return accelerations;
        }

        using tensor_vector_type = Eigen::Tensor<scalar_type, 3>;
        using tensor_scalar_type = Eigen::Tensor<scalar_type, 1>;

        // Tensor (grid) version
        tensor_scalar_type tensor_potential(
                const tensor_vector_type &x_grid,
                const tensor_vector_type &y_grid,
                const tensor_vector_type &z_grid) {
            std::array<tensor_vector_type, 3> grid{x_grid, y_grid, z_grid};
            tensor_scalar_type                potential(x_grid.dimension(0), y_grid.dimension(1), z_grid.dimension(2));

            GridCalculator<3, 1>::calculate(
                    *this,
                    grid,
                    [](MockModel &model, const auto &var_set) -> scalar_type {
                        return model.potential(var_set);
                    },
                    potential);
            return potential;
        }

        // Tensor (grid) version
        tensor_vector_type tensor_acceleration(
                const tensor_vector_type &x_grid,
                const tensor_vector_type &y_grid,
                const tensor_vector_type &z_grid) {
            std::array<tensor_vector_type, 3> grid{x_grid, y_grid, z_grid};
            tensor_vector_type                acceleration(x_grid.dimension(0), y_grid.dimension(1), z_grid.dimension(2));

            GridCalculator<3, 3>::calculate(
                    *this,
                    grid,
                    [](MockModel &model, const auto &var_set) -> vector_type {
                        return model.acceleration(var_set);
                    },
                    acceleration);
            return acceleration;
        }
    };

    using MockIMatrix = Eigen::Matrix<MockModel::scalar_type, Eigen::Dynamic, MockModel::Dim>;
    using MockOMatrix = Eigen::Matrix<MockModel::scalar_type, Eigen::Dynamic, 1>;

    class SeriesCalculatorTest : public ::testing::Test {
    protected:
        MockModel   model;
        MockIMatrix positions;

        SeriesCalculatorTest() : model(2), positions(MockIMatrix::Random(10, MockModel::Dim)) {}
    };

    TEST_F(SeriesCalculatorTest, SeriesPotentialTest) {
        MockOMatrix expected(positions.rows());
        for (int i = 0; i < positions.rows(); ++i) {
            expected(i, 0) = model.potential(positions.row(i));
        }

        MockOMatrix actual = model.potential_series(positions);
        ASSERT_TRUE(expected.isApprox(actual));
    }

    TEST_F(SeriesCalculatorTest, SeriesAccelerationTest) {
        MockIMatrix expected(positions.rows(), MockModel::Dim);
        for (int i = 0; i < positions.rows(); ++i) {
            expected.row(i) = model.acceleration(positions.row(i));
        }

        MockIMatrix actual = model.acceleration_series(positions);
        ASSERT_TRUE(expected.isApprox(actual));
    }

    TEST_F(SeriesCalculatorTest, SeriesPotentialEmptyTest) {
        MockIMatrix empty_positions(0, MockModel::Dim);
        MockOMatrix empty_actual = model.potential_series(empty_positions);
        ASSERT_EQ(empty_actual.size(), 0);
    }

    TEST_F(SeriesCalculatorTest, SeriesPotentialSinglePositionTest) {
        MockIMatrix single_position = MockIMatrix::Random(1, MockModel::Dim);
        MockOMatrix expected(1);
        expected(0, 0)     = model.potential(single_position.row(0));
        MockOMatrix actual = model.potential_series(single_position);
        ASSERT_TRUE(expected.isApprox(actual));
    }

    TEST_F(SeriesCalculatorTest, SeriesPotentialLargeInputTest) {
        MockIMatrix large_positions = MockIMatrix::Random(10000, MockModel::Dim);
        MockOMatrix expected(10000);
        for (int i = 0; i < large_positions.rows(); ++i) {
            expected(i, 0) = model.potential(large_positions.row(i));
        }
        MockOMatrix actual = model.potential_series(large_positions);
        ASSERT_TRUE(expected.isApprox(actual));
    }

    TEST_F(SeriesCalculatorTest, SeriesPotentialPositionsWithZerosTest) {
        MockIMatrix positions_with_zeros = MockIMatrix::Zero(10, MockModel::Dim);
        MockOMatrix expected(10);
        for (int i = 0; i < positions_with_zeros.rows(); ++i) {
            expected(i, 0) = model.potential(positions_with_zeros.row(i));
        }
        MockOMatrix actual = model.potential_series(positions_with_zeros);
        ASSERT_TRUE(expected.isApprox(actual));
    }

    class GridCalculatorTest : public ::testing::Test {
    protected:
        MockModel                     model;
        MockModel::tensor_vector_type x_grid, y_grid, z_grid;

        GridCalculatorTest() : model(2) {
            x_grid = MockModel::tensor_vector_type(10, 10, 10).setRandom();
            y_grid = MockModel::tensor_vector_type(10, 10, 10).setRandom();
            z_grid = MockModel::tensor_vector_type(10, 10, 10).setRandom();
        }
    };

    TEST_F(GridCalculatorTest, TensorPotentialTest) {
        MockModel::tensor_scalar_type expected(10, 10, 10);
        expected.setZero();
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                for (int k = 0; k < 10; ++k) {
                    Eigen::Vector3i index(i, j, k);
                    expected(index) = model.potential(Eigen::Vector3d{x_grid(index), y_grid(index), z_grid(index)});
                }
            }
        }

        MockModel::tensor_scalar_type actual = model.tensor_potential(x_grid, y_grid, z_grid);
        ASSERT_TRUE((actual - expected).abs().maxCoeff() < 1e-5);
    }

    TEST_F(GridCalculatorTest, TensorAccelerationTest) {
        MockModel::tensor_vector_type expected(10, 10, 10);
        expected.setZero();
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                for (int k = 0; k < 10; ++k) {
                    Eigen::Vector3i index(i, j, k);
                    expected(index) = model.acceleration(Eigen::Vector3d{x_grid(index), y_grid(index), z_grid(index)});
                }
            }
        }

        MockModel::tensor_vector_type actual = model.tensor_acceleration(x_grid, y_grid, z_grid);
        ASSERT_TRUE((actual - expected).abs().maxCoeff() < 1e-5);
    }

}// namespace
