load("@rules_cc//cc:defs.bzl", "cc_test")

alias(
    name = "gtest_main",
    actual = "@com_google_googletest//:gtest_main",
)

cc_test(
    name = "rotation_matrix_generator",
    srcs = ["test_rotation_matrix_generator.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:eigen",
        "@//odin:odin",
    ],
)

cc_test(
    name = "history",
    srcs = ["test_history.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:eigen",
        "@//odin:odin",
    ],
)

cc_test(
    name = "tracker",
    srcs = ["test_tracker.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:eigen",
        "@//odin:odin",
    ],
)

cc_test(
    name = "container_policy",
    srcs = ["test_container_policy.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:cereal",
        "@//odin:eigen",
        "@//odin:odin",
    ],
)

cc_test(
    name = "serialization_policy",
    srcs = ["test_serialization_policy.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:cereal",
        "@//odin:eigen",
        "@//odin:odin",
    ],
)

cc_test(
    name = "runge_kutta",
    srcs = ["test_runge_kutta.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:odin",
    ],
)

cc_test(
    name = "runge_kutta_embedded",
    srcs = ["test_runge_kutta_embedded.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:odin",
    ],
)

test_suite(
    name = "all_core",
    tests = [
        ":container_policy",
        ":history",
        ":rotation_matrix_generator",
        ":runge_kutta",
        ":tracker",
    ],
)
