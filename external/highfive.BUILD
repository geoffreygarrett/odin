load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

filegroup(
    name = "highfive_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:private"],
    # Private is default, but just to be explicit. Sources would only be
    # part of the cmake target, but not part of the cc_library target.
)

filegroup(
    name = "highfive_headers",
    srcs = glob(["include/**"]),
    visibility = ["//visibility:public"],
    # We make this public, so we have the option of packaging it
    # with our output.
)

# Even though eigen is header only, libraries like pagmo are easier
# to build with a cmake target. So we create a cmake target that
# performs the installation of the headers and cmake configs.
cmake(
    name = "highfive_cmake",
    lib_source = ":highfive_sources",
    out_include_dir = "include",
    out_headers_only = True,  # Flag variable to indicate that the library produces only headers
    install = True,
    deps = [
        "@com_github_hdfgroup_hdf5//:hdf5",
        "@com_github_eigen_eigen//:eigen",
    ],
    cache_entries = {
        "CMAKE_CXX_STANDARD": "20",
        "HIGHFIVE_USE_EIGEN": "ON",
        "HIGHFIVE_USE_BOOST": "OFF",
        "HIGHFIVE_PARALLEL_HDF5": "OFF",  # need to configure later in hdf5
        "HIGHFIVE_USE_XTENSOR": "OFF",
        "HIGHFIVE_USE_OPENCV": "OFF",
        "HIGHFIVE_USE_HALF_FLOAT": "OFF",
        "Eigen3_DIR": "$$EXT_BUILD_DEPS$$/eigen/share/eigen3/cmake",
        #        "Eigen3_DIR": "@com_github_eigen_eigen//:eigen/share/eigen3/cmake",
        # TODO: Integrate this later with threading options etc.
    },
    visibility = ["//visibility:public"],
    linkopts = ["-lz"],
)

cc_library(
    name = "highfive",
    hdrs = [":highfive_cmake"],
    visibility = ["//visibility:public"],
    includes = ["include"],
)
