load("@rules_cc//cc:defs.bzl", "cc_library")

#cc_library(
#    name = "autodiff",
#    srcs = [],
#    hdrs = glob(["autodiff/**/*.hpp"]),
#    includes = ['.'],
#    visibility = ["//visibility:public"],
#)
#cc_library(
#    name = "autodiff",
#    hdrs = glob(["autodiff/**/*.hpp"]),
#    includes = ["."],  # Include the parent directory
#    visibility = ["//visibility:public"],
#)

#cc_library(
#    name = "common",
#    hdrs = glob(["common/*.hpp"]),
#    deps = ["@com_github_eigen_eigen//:eigen"],
#    visibility = ["//visibility:public"],
#)
#
#cc_library(
#    name = "reverse",
#    hdrs = glob(["var.hpp", "reverse/**/*.hpp"]),
#    deps = ["@com_github_eigen_eigen//:eigen",
#            ":common"],
#    visibility = ["//visibility:public"],
#)
#
#cc_library(
#    name = "forward",
#    hdrs = glob(["dual.hpp", "real.hpp", "forward/**/*.hpp"]),
#    deps = [
#        "@com_github_eigen_eigen//:eigen",
#        ":common",
#    ],
#    visibility = ["//visibility:public"],
#)

cc_library(
    name = "autodiff",
    deps = ["//autodiff:reverse",
            "//autodiff:forward"],
    visibility = ["//visibility:public"],
    includes = ["."],
    hdrs = glob(["autodiff/**"]),
)

