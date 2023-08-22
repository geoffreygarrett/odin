load("@rules_cc//cc:defs.bzl", "cc_test")

alias(
    name = "gtest_main",
    actual = "@com_google_googletest//:gtest_main",
)

cc_test(
    name = "astrodynamics",
    srcs = ["test_astrodynamics.cpp"],
    deps = [
        ":gtest_main",
        "@//odin:odin",
    ],
)

test_suite(
    name = "all_domain",
    tests = [
        ":astrodynamics",
    ],
)
