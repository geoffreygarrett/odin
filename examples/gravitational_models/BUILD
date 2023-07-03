cc_binary(
    name = "tri_axial_ellipsoidal_1",
    srcs = ["tri_axial_ellipsoidal_1.cpp"],
    copts = ["-std=c++20"],
    defines = [
        "ODIN_USE_GLOG",
        "GLOG_CUSTOM_PREFIX_SUPPORT",
        "HIGHFIVE_USE_EIGEN=1",
        "H5_USE_EIGEN=1",  # -.-
    ],
    includes = ["$$EXT_BUILD_DEPS$$/highfive_cmake/include"],
    deps = [
        "//:gsl",
        "//:highfive",
        "//:odin",
        "//:tbb",
        "@com_github_google_glog//:glog",
        "@com_google_googletest//:gtest",
    ],
)

cc_binary(
    name = "tri_axial_ellipsoidal_2",
    srcs = ["tri_axial_ellipsoidal_2.cpp"],
    copts = ["-std=c++20"],
    defines = [
        "ODIN_USE_GLOG",
        "GLOG_CUSTOM_PREFIX_SUPPORT",
        "HIGHFIVE_USE_EIGEN=1",
        "H5_USE_EIGEN=1",  # -.-
    ],
    includes = ["$$EXT_BUILD_DEPS$$/highfive_cmake/include"],
    deps = [
        "//:gsl",
        "//:highfive",
        "//:odin",
        "//:tbb",
        #        "//:shtns",
        "@com_github_google_glog//:glog",
        #        "@com_google_googletest//:gtest",
    ],
)

cc_binary(
    name = "tri_axial_ellipsoidal_orbit",
    srcs = ["tri_axial_ellipsoidal_orbit.cpp"],
    defines = [
        "ODIN_USE_GLOG",
        "GLOG_CUSTOM_PREFIX_SUPPORT",
        "HIGHFIVE_USE_EIGEN=1",
        "H5_USE_EIGEN=1",  # -.-
    ],
    includes = ["$$EXT_BUILD_DEPS$$/highfive_cmake/include"],
    deps = [
        "//:gsl",
        "//:highfive",
        "//:odin",
        "//:tbb",
        #        "//:shtns",
        "@com_github_google_glog//:glog",
        #        "@com_google_googletest//:gtest",
    ],
)

cc_binary(
    name = "spherical_harmonics",
    srcs = ["spherical_harmonics.cpp"],
    defines = [
        "ODIN_USE_GLOG",
        "GLOG_CUSTOM_PREFIX_SUPPORT",
        "HIGHFIVE_USE_EIGEN=1",
        "H5_USE_EIGEN=1",  # -.-
    ],
    includes = ["$$EXT_BUILD_DEPS$$/highfive_cmake/include"],
    deps = [
        "//:gsl",
        "//:highfive",
        "//:odin",
        "//:tbb",
        #        "//:shtns",
        "@com_github_google_glog//:glog",
        #        "@com_google_googletest//:gtest",
    ],
)

cc_binary(
    name = "sh_jgm_3",
    srcs = ["sh_jgm_3.cpp"],
    defines = [
        "ODIN_USE_GLOG",
        "GLOG_CUSTOM_PREFIX_SUPPORT",
        "HIGHFIVE_USE_EIGEN=1",
        "H5_USE_EIGEN=1",  # -.-
    ],
    includes = ["$$EXT_BUILD_DEPS$$/highfive_cmake/include"],
    deps = [
        "//:gsl",
        "//:highfive",
        "//:odin",
        "//:tbb",
        #        "//:shtns",
        "@com_github_google_glog//:glog",
        #        "@com_google_googletest//:gtest",
    ],
)

genrule(
    name = "run_tri_axial_ellipsoidal_2",
    outs = [
        "point_mass.h5",
        "ellipsoid.h5",
    ],  # or whatever file your cc_binary outputs
    cmd = "$(location :tri_axial_ellipsoidal_2)",
    tools = [":tri_axial_ellipsoidal_2"],
)

py_binary(
    name = "plot",
    srcs = ["plot.py"],
    data = [
        ":run_tri_axial_ellipsoidal_2",  # your plot.py script can then use output_file
    ],
    deps = [],
)

#py_binary(
#    name = "plot",
#    srcs = ["plot.py"],
#    data = [
#        ":tri_axial_ellipsoidal_2",
#    ],
#    deps = [],
#)
