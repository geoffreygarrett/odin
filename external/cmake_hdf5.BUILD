package(default_visibility = ["//visibility:public"])

# other example https://github.com/ibab/laminate/blob/master/bazel/hdf5.BUILD
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

filegroup(
    name = "hdf5_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:private"],
    # Private is default, but just to be explicit. Sources would only be
    # part of the cmake target, but not part of the cc_library target.
)

filegroup(
    name = "hdf5_headers",
    srcs = glob(["include/**"]),
    visibility = ["//visibility:public"],
    # We make this public, so we have the option of packaging it
    # with our output.
)

cmake(
    name = "hdf5_cmake",
    lib_source = ":hdf5_sources",
    out_include_dir = "include",
    install = True,
    # Bazel focuses on hermetic and to have build reproducible. This is why all items such as
    # __DATE__, __TIMESTAMP__, and __TIME__ are replaced with "redacted" to avoid the timestamp
    # changing for each build, the only problem is if these are used In C code, it will not compile.
    # We need to escape the quotes correctly as seen here.
    cache_entries = {
        "CMAKE_CXX_STANDARD": "20",
        #        "HIGHFIVE_USE_EIGEN": "ON",
        #        "HIGHFIVE_USE_BOOST": "OFF",
        #        "HIGHFIVE_PARALLEL_HDF5": "ON",
        #        "HIGHFIVE_USE_XTENSOR": "OFF",
        #        "HIGHFIVE_USE_OPENCV": "OFF",
        #        "HIGHFIVE_USE_HALF_FLOAT": "OFF",
        # TODO: Integrate this later with threading options etc.
    },
    visibility = ["//visibility:public"],
)

cc_library(
    name = "hdf5",
    deps = [":hdf5_cmake"],
    visibility = ["//visibility:public"],
    includes = ["include"],
)
