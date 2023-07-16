#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <odin/vectorize.hpp>
#include <unsupported/Eigen/CXX11/Tensor>
#include <utility>

//numpy.linspace
//        numpy.linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None, axis=0)[source]
//        Return evenly spaced numbers over a specified interval.
//
//        Returns num evenly spaced samples, calculated over the interval [start, stop].
//
//                The endpoint of the interval can optionally be excluded.
//
//                Changed in version 1.16.0: Non-scalar start and stop are now supported.
//
//                  Changed in version 1.20.0: Values are rounded towards -inf instead of 0 when an integer dtype is specified. The old behavior can still be obtained with np.linspace(start, stop, num).astype(int)
//
//                          Parameters:
//    startarray_like
//    The starting value of the sequence.
//
//    stoparray_like
//    The end value of the sequence, unless endpoint is set to False. In that case, the sequence consists of all but the last of num + 1 evenly spaced samples, so that stop is excluded. Note that the step size changes when endpoint is False.
//
//                                                                                                                                  numint, optional
//                Number of samples to generate. Default is 50. Must be non-negative.
//
//                  endpointbool, optional
//        If True, stop is the last sample. Otherwise, it is not included. Default is True.
//
//                                           retstepbool, optional
//        If True, return (samples, step), where step is the spacing between samples.
//
//                                dtypedtype, optional
//        The type of the output array. If dtype is not given, the data type is inferred from start and stop. The inferred dtype will never be an integer; float is chosen even if the arguments would produce an array of integers.
//
//        New in version 1.9.0.
//
//        axisint, optional
//                The axis in the result to store the samples. Relevant only if start or stop are array-like. By default (0), the samples will be along a new axis inserted at the beginning. Use -1 to get an axis at the end.
//
//                                                                             New in version 1.16.0.
//
//                                                                             Returns:
//samplesndarray
//There are num equally spaced samples in the closed interval [start, stop] or the half-open interval [start, stop) (depending on whether endpoint is True or False).
//
//stepfloat, optional
//Only returned if retstep is True
//
//Size of spacing between samples.
//
//See also
//
//arange
//Similar to linspace, but uses a step size (instead of the number of samples).
//
//geomspace
//Similar to linspace, but with numbers spaced evenly on a log scale (a geometric progression).
//
//logspace
//Similar to geomspace, but with the end points specified as logarithms.
//
//How to create arrays with regularly-spaced values
//Examples
//
//np.linspace(2.0, 3.0, num=5)
//array([2.  , 2.25, 2.5 , 2.75, 3.  ])
//np.linspace(2.0, 3.0, num=5, endpoint=False)
//array([2. ,  2.2,  2.4,  2.6,  2.8])
//np.linspace(2.0, 3.0, num=5, retstep=True)
//(array([2.  ,  2.25,  2.5 ,  2.75,  3.  ]), 0.25)
using namespace odin;


class LinspaceTest : public ::testing::Test {
protected:
    // Create lambda functions that capture the correct overloads
    std::function<Eigen::VectorXd(double, double, int, bool, bool)> linspace =
            [](double start, double end, int num_points, bool include_endpoint, bool parallel) -> Eigen::VectorXd {
        return eigen::linspace(start, end, num_points, include_endpoint, parallel);
    };

    void run_linspace_tests(double start, double end, int num_points, bool include_endpoint, bool parallel,
                            const Eigen::VectorXd &expected_values) {
        Eigen::VectorXd vec = linspace(start, end, num_points, include_endpoint, parallel);
        ASSERT_EQ(vec.size(), expected_values.size());
        for (Eigen::Index i = 0; i < vec.size(); i++) {
            ASSERT_NEAR(vec(i), expected_values(i), 1e-9);// Adjust tolerance as necessary
        }
    }
};

TEST_F(LinspaceTest, BasicTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(5, 1., 5.);
    run_linspace_tests(1., 5., 5, true, false, expected);
}

TEST_F(LinspaceTest, NegativeStartTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(5, -2., 2.);
    run_linspace_tests(-2., 2., 5, true, false, expected);
}

TEST_F(LinspaceTest, ExcludeEndpointTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(5, -2., 2. - 4. / 5.);
    run_linspace_tests(-2., 2., 5, false, false, expected);
}

TEST_F(LinspaceTest, LargeNumTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(1e6 + 1, 0., 1e6);
    run_linspace_tests(0., 1e6, 1e6 + 1, true, true, expected);// Run in parallel
}

TEST_F(LinspaceTest, OnePointTest) {
    Eigen::VectorXd expected(1);
    expected << 1.0;
    run_linspace_tests(1., 5., 1, true, false, expected);
}

TEST_F(LinspaceTest, ZeroPointsTest) {
    Eigen::VectorXd expected(0);
    run_linspace_tests(1., 5., 0, true, false, expected);
}

TEST_F(LinspaceTest, PositiveAndNegativeValuesTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(5, -3., 2.);
    run_linspace_tests(-3., 2., 5, true, false, expected);
}

TEST_F(LinspaceTest, VeryLargeValuesTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(5, 1e12, 1e13);
    run_linspace_tests(1e12, 1e13, 5, true, true, expected);// Run in parallel
}

TEST_F(LinspaceTest, VerySmallValuesTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(5, 1e-12, 1e-11);
    run_linspace_tests(1e-12, 1e-11, 5, true, false, expected);
}

TEST_F(LinspaceTest, EndpointAtZeroTest) {
    Eigen::VectorXd expected = Eigen::VectorXd::LinSpaced(5, -2., 0.);
    run_linspace_tests(-2., 0., 5, true, false, expected);
}

class LogspaceTest : public ::testing::Test {
protected:
    // Create lambda functions that capture the correct overloads
    std::function<Eigen::VectorXd(double, double, int, bool, bool)> logspace =
            [](double start, double end, int num_points, bool include_endpoint, bool parallel) -> Eigen::VectorXd {
        return eigen::logspace(start, end, num_points, include_endpoint, parallel);
    };

    static Eigen::VectorXd expected_logspace(double start, double end, int num_points, bool endpoint, bool parallel = false) {
        Eigen::VectorXd output(num_points);
        if (num_points > 1) {
            double exp_start = std::pow(10, start);
            double exp_end   = std::pow(10, end);
            double ratio     = endpoint ? std::pow(exp_end / exp_start, 1.0 / (num_points - 1)) : std::pow(exp_end / exp_start, 1.0 / num_points);
            if (parallel) {
                //#pragma omp parallel for
                for (int i = 0; i < num_points; ++i) {
                    output(i) = exp_start * std::pow(ratio, i);
                }
            } else {
                for (int i = 0; i < num_points; ++i) {
                    output(i) = exp_start * std::pow(ratio, i);
                }
            }
            if (!endpoint) {
                output(num_points - 1) = std::pow(10, end);
            }
        } else if (num_points == 1) {
            output(0) = std::pow(10, start);
        }
        return output;
    }

    void run_logspace_tests(double start, double end, int num_points, bool include_endpoint, bool parallel,
                            const Eigen::VectorXd &expected_values) {
        Eigen::VectorXd vec = logspace(start, end, num_points, include_endpoint, parallel);
        ASSERT_EQ(vec.size(), expected_values.size());
        for (Eigen::Index i = 0; i < vec.size(); i++) {
            if (std::isinf(vec(i)) && std::isinf(expected_values(i))) {
                // Both values are inf, hence they are considered equal.
                continue;
            } else if (std::isinf(vec(i)) || std::isinf(expected_values(i))) {
                FAIL() << "One of the values is inf, but the other is not.";
            } else if (expected_values(i) == 0.0) {
                ASSERT_NEAR(vec(i), 0.0, 1e-9);
            } else {
                ASSERT_NEAR(vec(i), expected_values(i), 1e-9);
            }
        }
    }
};

TEST_F(LogspaceTest, BasicTest) {
    Eigen::VectorXd expected = expected_logspace(1., 5., 5, true, false);
    run_logspace_tests(1., 5., 5, true, false, expected);
}

TEST_F(LogspaceTest, LargeNumTest) {
    Eigen::VectorXd expected = expected_logspace(0., 1e6, 1e6 + 1, true, true);
    run_logspace_tests(0., 1e6, 1e6 + 1, true, true, expected);
}

TEST_F(LogspaceTest, OnePointTest) {
    Eigen::VectorXd expected(1);
    expected << 1.0;// start is 1, meaning e^1
    run_logspace_tests(1., 5., 1, true, false, expected);
}

TEST_F(LogspaceTest, ZeroPointsTest) {
    Eigen::VectorXd expected(0);
    run_logspace_tests(1., 5., 0, true, false, expected);
}

TEST_F(LogspaceTest, VeryLargeValuesTest) {
    Eigen::VectorXd expected = expected_logspace(12, 13, 5, true, true);
    run_logspace_tests(12, 13, 5, true, true, expected);
}

TEST_F(LogspaceTest, VerySmallValuesTest) {
    Eigen::VectorXd expected = expected_logspace(-12, -11, 5, true, false);
    run_logspace_tests(-12, -11, 5, true, false, expected);
}

TEST_F(LogspaceTest, EndpointAtOneTest) {
    Eigen::VectorXd expected = expected_logspace(-2., 0., 5, true, false);
    run_logspace_tests(-2., 0., 5, true, false, expected);
}

class GeomspaceTest : public ::testing::Test {
protected:
    std::function<Eigen::VectorXd(double, double, int, bool, bool)> geomspace =
            [](double start, double end, int num_points, bool include_endpoint, bool parallel) -> Eigen::VectorXd {
        return eigen::geomspace(start, end, num_points, include_endpoint, parallel);
    };

    static Eigen::VectorXd expected_geomspace(double start, double end, int num_points, bool endpoint, bool parallel = false) {
        Eigen::VectorXd output(num_points);
        if (num_points > 1) {
            double ratio = endpoint ? std::pow(end / start, 1.0 / (num_points - 1)) : std::pow(end / start, 1.0 / num_points);
            if (parallel) {
                //#pragma omp parallel for
                for (int i = 0; i < num_points; ++i) {
                    output(i) = start * std::pow(ratio, i);
                }
            } else {
                for (int i = 0; i < num_points; ++i) {
                    output(i) = start * std::pow(ratio, i);
                }
            }
            if (!endpoint) {
                output(num_points - 1) = end;
            }
        } else if (num_points == 1) {
            output(0) = start;
        }
        return output;
    }

    void run_geomspace_tests(double start, double end, int num_points, bool include_endpoint, bool parallel,
                            const Eigen::VectorXd &expected_values) {
        Eigen::VectorXd vec = geomspace(start, end, num_points, include_endpoint, parallel);
        ASSERT_EQ(vec.size(), expected_values.size());
        for (Eigen::Index i = 0; i < vec.size(); i++) {
            if (std::isinf(vec(i)) && std::isinf(expected_values(i))) {
                continue;
            } else if (std::isinf(vec(i)) || std::isinf(expected_values(i))) {
                FAIL() << "One of the values is inf, but the other is not.";
            } else if (std::isnan(vec(i)) && std::isnan(expected_values(i))) {
                continue;
            } else if (std::isnan(vec(i)) || std::isnan(expected_values(i))) {
                FAIL() << "One of the values is NaN, but the other is not.";
            } else if (expected_values(i) == 0.0) {
                ASSERT_NEAR(vec(i), 0.0, 1e-9);
            } else {
                ASSERT_NEAR(vec(i), expected_values(i), 1e-9);
            }
        }
    }
};

TEST_F(GeomspaceTest, BasicTest) {
    Eigen::VectorXd expected = expected_geomspace(1., 5., 5, true, false);
    run_geomspace_tests(1., 5., 5, true, false, expected);
}

TEST_F(GeomspaceTest, LargeNumTest) {
    Eigen::VectorXd expected = expected_geomspace(0., 1e6, 1e6 + 1, true, true);
    run_geomspace_tests(0., 1e6, 1e6 + 1, true, true, expected);
}

TEST_F(GeomspaceTest, OnePointTest) {
    Eigen::VectorXd expected(1);
    expected << 1.0;
    run_geomspace_tests(1., 5., 1, true, false, expected);
}

TEST_F(GeomspaceTest, ZeroPointsTest) {
    Eigen::VectorXd expected(0);
    run_geomspace_tests(1., 5., 0, true, false, expected);
}

TEST_F(GeomspaceTest, VeryLargeValuesTest) {
    Eigen::VectorXd expected = expected_geomspace(12, 13, 5, true, true);
    run_geomspace_tests(12, 13, 5, true, true, expected);
}

TEST_F(GeomspaceTest, VerySmallValuesTest) {
    Eigen::VectorXd expected = expected_geomspace(-12, -11, 5, true, false);
    run_geomspace_tests(-12, -11, 5, true, false, expected);
}

TEST_F(GeomspaceTest, EndpointAtOneTest) {
    Eigen::VectorXd expected = expected_geomspace(-2., 0., 5, true, false);
    run_geomspace_tests(-2., 0., 5, true, false, expected);
}

class MeshgridTest : public testing::Test {
protected:
    template<typename T>
    std::tuple<Eigen::VectorXd, Eigen::VectorXd, std::array<Eigen::Tensor<double, 2>, 2>, std::array<Eigen::Tensor<double, 2>, 2>>
    create_test_inputs_and_outputs(T start_x, T end_x, T start_y, T end_y, int num_points) {
        Eigen::VectorXd x_vec = eigen::linspace(start_x, end_x, num_points);
        Eigen::VectorXd y_vec = eigen::linspace(start_y, end_y, num_points);

        // Converting VectorXd to Tensor for assignment
        Eigen::Tensor<double, 1> x = Eigen::TensorMap<Eigen::Tensor<const double, 1>>(x_vec.data(), x_vec.size());
        Eigen::Tensor<double, 1> y = Eigen::TensorMap<Eigen::Tensor<const double, 1>>(y_vec.data(), y_vec.size());

        Eigen::Tensor<double, 2> expected_x(y.size(), x.size());
        for (int i = 0; i < y.size(); ++i) {
            expected_x.chip(i, 0) = x;
        }

        Eigen::Tensor<double, 2> expected_y(y.size(), x.size());
        for (int i = 0; i < x.size(); ++i) {
            expected_y.chip(i, 1) = y;
        }
        std::array<Eigen::Tensor<double, 2>, 2> expected_xy = {expected_x, expected_y};

        Eigen::array<int, 2>                    transpose_dims = {1, 0};
        std::array<Eigen::Tensor<double, 2>, 2> expected_ij    = {
                expected_x.shuffle(transpose_dims),
                expected_y.shuffle(transpose_dims)};

        return std::make_tuple(x_vec, y_vec, expected_xy, expected_ij);
    }

    template<typename T>
    void run_meshgrid_tests(T start_x, T end_x, T start_y, T end_y, int num_points) {
        auto [x, y, expected_xy, expected_ij] = create_test_inputs_and_outputs(start_x, end_x, start_y, end_y, num_points);
        auto [X_xy, Y_xy]                     = meshgrid(Indexing::xy, true, x, y);
        auto [X_ij, Y_ij]                     = meshgrid(Indexing::ij, true, x, y);

        EXPECT_TRUE(are_tensors_equal(X_xy, expected_xy[0]));
        EXPECT_TRUE(are_tensors_equal(Y_xy, expected_xy[1]));
        EXPECT_TRUE(are_tensors_equal(X_ij, expected_ij[0]));
        EXPECT_TRUE(are_tensors_equal(Y_ij, expected_ij[1]));
    }

    static bool are_tensors_equal(const Eigen::Tensor<double, 2> &tensor1, const Eigen::Tensor<double, 2> &tensor2) {
        Eigen::Tensor<bool, 0> tensorResult = ((tensor1 - tensor2).abs() >= 1e-5).any();
        return !tensorResult();
    }
};

TEST_F(MeshgridTest, BasicTest) {
    run_meshgrid_tests(1., 3., 4., 6., 3);
}

TEST_F(MeshgridTest, LargeNumTest) {
    run_meshgrid_tests(0., 1e3, 0., 1e3, 1e3 + 1);
}

TEST_F(MeshgridTest, NonEqualVectorSizes) {
    run_meshgrid_tests(-1., 1., 1., 3., 4);
}

TEST_F(MeshgridTest, NegativeRangeInputs) {
    run_meshgrid_tests(-2., 2., -3., 3., 7);
}

TEST_F(MeshgridTest, FloatingPointValues) {
    run_meshgrid_tests(0.5, 1.5, 2.5, 3.5, 3);
}

TEST_F(MeshgridTest, LargeNumbers) {
    run_meshgrid_tests(1e6, 1e6 + 5., 1e7, 1e7 + 5., 5);
}

TEST_F(MeshgridTest, NonUniformlySpacedInputs) {
    run_meshgrid_tests(1., 8., 10., 40., 4);
}

TEST_F(MeshgridTest, ZeroToOneInput) {
    run_meshgrid_tests(0., 1., 0., 1., 3);
}

TEST_F(MeshgridTest, NegativeToOneInput) {
    run_meshgrid_tests(-1., 1., -1., 1., 5);
}