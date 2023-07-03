cc_binary(
    name = "main",
    srcs = ["main.cpp"],
    defines = [
        "ODIN_USE_GLOG",
        "GLOG_CUSTOM_PREFIX_SUPPORT",
    ],
    deps = [
        "//:gsl",
        "//:highfive",
        "//:odin",
        "//:tbb",
        "@com_github_google_glog//:glog",
        "@com_google_googletest//:gtest",
    ],
)
