load("@rules_cc//cc:defs.bzl", "cc_test")

alias(
    name = "gtest_main",
    actual = "@com_google_googletest//:gtest_main",
)

cc_test(
    name = "autodiff",
    srcs = ["test_autodiff.cpp"],
    deps = [
        ":gtest_main",
        "//:autodiff",  # their pacakge already contains eigen as a dependency
    ],
)

cc_test(
    name = "cereal",
    srcs = ["test_cereal.cpp"],
    deps = [
        ":gtest_main",
        "//:cereal",
        "//:odin",
        "//:eigen",
        # TODO: This should be a test for policy/serialization
    ],
)

#cc_library(
#    name = "pagmo",
#    srcs = glob([
#        "miniconda/envs/odin-env/**/*/libpagmo*.so*",
#        "miniconda/envs/odin-env/**/*.so*",
#        "miniconda/envs/odin-env/**/*.a",
#    ]),
#    hdrs = glob([
#        "miniconda/envs/odin-env/include/pagmo/**/*.hpp",
#        "miniconda/envs/odin-env/include/pagmo/**/*.h",
#        "miniconda/envs/odin-env/include/**/*.h",
#        "miniconda/envs/odin-env/include/**/*.hpp",
#    ]),
#    #    deps = [":boost", ":tbb"],  # Pagmo depends on Boost
#    includes = ["miniconda/envs/odin-env/include"],
#    linkopts = ["-Wl,-rpath,miniconda/envs/odin-env/lib"],
#    visibility = ["//visibility:public"],
#    deps = [
#        "@conda_prefix//:include",
#        "@conda_prefix//:lib",
#    ],
#)

#cc_import(
#    name = "pagmo",
#    hdrs = glob([
#        "include/pagmo/**/*.hpp",
#        "include/pagmo/**/*.h",
#    ]),
#    data = [
#        "@conda//:prefix",
#    ],
#    visibility = ["//visibility:public"],
#)
#cc_library(
#    name = "pagmo",
#    hdrs = glob([
#        "include/pagmo/**/*.hpp",
#        "include/pagmo/**/*.h",
#    ]),
#    includes = ["include"],
#    visibility = ["//visibility:public"],
#    deps = [
#        "@conda//:prefix",
#    ],
#)
#cc_import(
#    name = "pagmo",
#    hdrs = glob([
#        "include/pagmo/**/*.hpp",
#        "include/pagmo/**/*.h",
#    ]),
#    data = [
#        "@conda//:prefix",
#    ],
#    #    shared_library = "lib/libpagmo.so",
#    visibility = ["//visibility:public"],
#)

#cc_test(
#    name = "test_pagmo",
#    srcs = ["test_pagmo.cpp"],
#    includes = ["include"],
#    deps = [
#        ":gtest_main",
#        "@//:pagmo",
#    ],
#)
##        "@odin_conda_environment//:pagmo",

#        "//:cereal",
#        "//:odin",
#        "//:eigen",
#        # TODO: This should be a test for policy/serialization
cc_test(
    name = "tbb",
    srcs = ["test_tbb.cpp"],
    deps = [
        ":gtest_main",
        "//:tbb",
    ],
)

test_suite(
    name = "all_external",
    tests = [
        ":autodiff",
        ":cereal",
        ":tbb",
        #        ":pagmo",
        #        ":gsl",
    ],
)
