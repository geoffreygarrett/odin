#include <Eigen/Dense>
#include <include/gtest/gtest.h>
#include <odin/core/frames/rotation_matrix_generator.hpp>

namespace {
    // Google Test framework allows us to define parameterized tests
    // Here we create a fixture class for such tests.
    class RotationMatrixTest : public testing::TestWithParam<std::tuple<Eigen::Vector2d, Eigen::Vector2d>> {
    };

    //    // Now we define a test for 2D case using the fixture
    //    TEST_P(RotationMatrixTest, Rotation2D) {
    //        Eigen::Vector2d p = std::get<0>(GetParam());
    //        Eigen::Vector2d expected = std::get<1>(GetParam());
    //
    //        auto R = make_rot_mat_lambda<axis::align<axis::x, 0>, axis::align<axis::x, 0>>()(p);
    //        Eigen::Vector2d result = R * p;
    //
    //        // Use Google Test's built-in support for Eigen to compare the resulting matrix with expected
    //        EXPECT_TRUE(result.isApprox(expected));
    //    }

    // Now define some parameters for the 2D test
    //    INSTANTIATE_TEST_SUITE_P(
    //            Rotation2DTest,
    //            RotationMatrixTest,
    //            testing::Values(
    //                    std::make_tuple(Eigen::Vector2d(1, 0), Eigen::Vector2d(1, 0)),
    //                    std::make_tuple(Eigen::Vector2d(0, 1), Eigen::Vector2d(0, -1)),
    //                    std::make_tuple(Eigen::Vector2d(-1, 0), Eigen::Vector2d(-1, 0)),
    //                    std::make_tuple(Eigen::Vector2d(0, -1), Eigen::Vector2d(0, 1))
    //            )
    //    );

    // Similarly, we can create a parameterized test for the 3D case
    class RotationMatrix3DTest
        : public testing::TestWithParam<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Matrix3d>> {
    };

    TEST_P(RotationMatrix3DTest, Rotation3D) {
        Eigen::Vector3d p = std::get<0>(GetParam());
        Eigen::Vector3d v = std::get<1>(GetParam());
        Eigen::Matrix3d expected = std::get<2>(GetParam());

        auto R = make_rot_mat_lambda<
                frame::align<frame::axis::x, 0>,
                frame::align<frame::axis::y, 1>>()(p, v);
        Eigen::Matrix3d result = R;

        EXPECT_TRUE(result.isApprox(expected));
    }

    Eigen::Matrix3d calculate_expected_rotation_matrix(const Eigen::Vector3d &p, const Eigen::Vector3d &v) {
        Eigen::Vector3d z = p.cross(v).normalized();
        Eigen::Matrix3d expected;
        expected.col(0) = p.normalized();
        expected.col(1) = v.normalized();
        expected.col(2) = z;
        return expected;
    }

    // TODO: this needs to cover every test matrix, but the setup to do so, is a little excessive right now.
    INSTANTIATE_TEST_SUITE_P(
            Rotation3DTest,
            RotationMatrix3DTest,
            testing::Values(
                    // Aligning primary vector with x-axis, secondary vector with y-axis
                    std::make_tuple(
                            Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0),
                            calculate_expected_rotation_matrix(Eigen::Vector3d(1, 0, 0),
                                                               Eigen::Vector3d(0, 1, 0))),
                    // Aligning primary vector with x-axis, secondary vector with z-axis
                    std::make_tuple(
                            Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 0, 1),
                            calculate_expected_rotation_matrix(Eigen::Vector3d(1, 0, 0),
                                                               Eigen::Vector3d(0, 0, 1))),

                    // Aligning primary vector with y-axis, secondary vector with x-axis
                    std::make_tuple(
                            Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(1, 0, 0),
                            calculate_expected_rotation_matrix(Eigen::Vector3d(0, 1, 0),
                                                               Eigen::Vector3d(1, 0, 0))),
                    // Aligning primary vector with y-axis, secondary vector with z-axis
                    std::make_tuple(
                            Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(0, 0, 1),
                            calculate_expected_rotation_matrix(Eigen::Vector3d(0, 1, 0),
                                                               Eigen::Vector3d(0, 0, 1))),

                    // Aligning primary vector with z-axis, secondary vector with x-axis
                    std::make_tuple(
                            Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(1, 0, 0),
                            calculate_expected_rotation_matrix(Eigen::Vector3d(0, 0, 1),
                                                               Eigen::Vector3d(1, 0, 0))),
                    // Aligning primary vector with z-axis, secondary vector with y-axis
                    std::make_tuple(
                            Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(0, 1, 0),
                            calculate_expected_rotation_matrix(Eigen::Vector3d(0, 0, 1),
                                                               Eigen::Vector3d(0, 1, 0))),

                    // Some additional test cases
                    std::make_tuple(
                            Eigen::Vector3d(1, 1, 0).normalized(), Eigen::Vector3d(0, 0, 1),
                            calculate_expected_rotation_matrix(Eigen::Vector3d(1, 1, 0).normalized(),
                                                               Eigen::Vector3d(0, 0, 1)))));
}// namespace
