load("@rules_cc//cc:defs.bzl", "cc_test")

alias(
    name = "gtest_main",
    actual = "@com_google_googletest//:gtest_main",
)

#cc_test(
#    name = "tbb",
#    srcs = ["test_tbb.cpp"],
#    deps = [
#        ":gtest_main",
#        "@//odin:tbb",
#    ],
#)

test_suite(
    name = "all_models",
    tests = [
        #        ":autodiff",
        #        ":cereal",
        #        ":tbb",
        #        ":gsl",
        #        ":pagmo2"
    ],
)
