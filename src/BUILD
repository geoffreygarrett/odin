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
