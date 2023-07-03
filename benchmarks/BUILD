load("@rules_cc//cc:defs.bzl", "cc_test")

cc_binary(
    name = "matrix_inversion",
    srcs = ["benchmark_matrix_inversion.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        "//:eigen",
        "//:odin",
        "@com_google_benchmark//:benchmark_main",
    ],
)
