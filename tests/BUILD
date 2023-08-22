load("@rules_cc//cc:defs.bzl", "cc_test")

alias(
    name = "gtest_main",
    actual = "@com_google_googletest//:gtest_main",
)

cc_test(
    name = "test_hello",
    srcs = ["test_hello.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_autodiff",
    srcs = ["test_autodiff.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:autodiff",
    ],
)

cc_test(
    name = "test_container_policy",
    srcs = ["test_container_policy.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_doppler_shift_measurement",
    srcs = ["test_doppler_shift_measurement.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_eigen",
    srcs = ["test_eigen.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_filters",
    srcs = ["test_filters.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_history",
    srcs = ["test_history.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_inverse",
    srcs = ["test_inverse.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_kalman",
    srcs = ["test_kalman.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_measurement",
    srcs = ["test_measurement.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_rotation_matrix_lambda",
    srcs = ["test_rotation_matrix_lambda.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_runge_kutta",
    srcs = ["test_runge_kutta.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_serialization_policy",
    srcs = ["test_serialization_policy.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_sigma_points",
    srcs = ["test_sigma_points.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_state",
    srcs = ["test_state.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_stat_handler",
    srcs = ["test_state_handler.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_tbb",
    srcs = ["test_tbb.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_tracker",
    srcs = ["test_tracker.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_tri_axial_ellipsoid_homogeneous_mass",
    srcs = ["test_tri_axial_ellipsoid_homogeneous_mass.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_ukf",
    srcs = ["test_ukf.cpp"],
    deps = [
        ":gtest_main",
    ],
)

cc_test(
    name = "test_tree",
    srcs = ["test_tree.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:odin",
    ],
)

cc_test(
    name = "test_vectorize",
    srcs = ["test_vectorize.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:odin",
    ],
)

cc_test(
    name = "test_parallel",
    srcs = ["test_parallel.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:odin",
    ],
)

cc_test(
    name = "test_linalg",
    srcs = ["test_linalg/test_basics.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:odin",
    ],
)

# All tests target
test_suite(
    name = "all_tests",
    tests = [
        "@//odin/tests/test_core:all_core",
        "@//odin/tests/test_external:all_external",
        "@//odin/tests/test_models:all_models",
        "@//odin/tests/test_domain:all_domain",
        #        ":test_autodiff",
        #        ":test_container_policy",
        #        ":test_doppler_shift_measurement",
        #        ":test_eigen",
        #        ":test_filters",
        #        ":test_hello",
        #        ":test_history",
        #        ":test_inverse",
        #        ":test_kalman",
        #        ":test_measurement",
        #        ":test_rotation_matrix_lambda",
        #        ":test_runge_kutta",
        #        ":test_serialization_policy",
        #        ":test_sigma_points",
        #        ":test_stat_handler",
        #        ":test_state",
        #        ":test_tbb",
        #        ":test_tracker",
        #        ":test_tri_axial_ellipsoid_homogeneous_mass",
        #        ":test_ukf",
    ],
)
