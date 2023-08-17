cc_binary(
    name = "main",
    srcs = ["main.cpp"],
    copts = select({
        "//:without_logging": [],
        "//conditions:default": ["-DODIN_USE_GLOG -DGLOG_CUSTOM_PREFIX_SUPPORT"],
    }),
    visibility = ["//visibility:public"],
    deps = select(
        {
            "//:without_logging": [],
            "//conditions:default": ["//:glog"],
        },
    ) + [
        "//:autodiff",
        "//:benchmark",
        "//:cereal",
        "//:eigen",
        "//:gsl",
        "//:odin",
        "//:tbb",
        #        "@odin_conda_environment//:pag",
    ],
)

load("@rules_cuda//cuda:defs.bzl", "cuda_library")

#cc_binary(
#    name = "polyhedral_gravity_test",
#    srcs = ["@odin//:polyhedral_gravity_main"],
#    deps = ["@com_github_esa_polyhedral_gravity_model//:polyhedral_gravity"],
#)

# thrust have code marked as __global__ or __kernel__
# so source for thrust should be built with cuda_library first
cuda_library(
    name = "thrust_cu",
    srcs = ["thrust.cu"],
    copts = ["-allow-unsupported-compiler"],
    defines = ["THRUST_IGNORE_CUB_VERSION_CHECK"],
    deps = ["@com_github_nvidia_thrust//:thrust"],
    #    deps = ["@local_cuda//:thrust"],
)

# Run with 'bazel run //examples/thrust:main'
#cc_binary(
#    name = "thrust_test",
#    deps = [":thrust_cu"],
#)

cc_binary(
    name = "main2",
    #    srcs = ["main.cu"],
    visibility = ["//visibility:public"],
    deps = [
        ":thrust_cu",
        "@local_cuda//:thrust",
        #        "//:odin",
        #        #        "@com_github_nvidia_thrust//:thrust",
        #        "@local_cuda//:thrust",
        #        #        "@com_github_esa_polyhedral_gravity_model//:polyhedral_gravity",
    ],
)
