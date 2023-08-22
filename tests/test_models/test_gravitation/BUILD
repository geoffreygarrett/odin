cc_binary(
    name = "tri_axial_ellipsoidal",
    srcs = ["test_tri_axial_ellipsoid.cpp"],
    defines = [
        "ODIN_USE_GLOG",
        "GLOG_CUSTOM_PREFIX_SUPPORT",
    ],
    deps = [
        "@//odin:gsl",
        "@//odin:odin",
        "@com_github_google_glog//:glog",
        "@com_google_googletest//:gtest",
    ],
)
