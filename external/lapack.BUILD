load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

cmake(
    name = "lapack_cmake",
    lib_source = ":lapack_sources",
    out_include_dir = "include",
    out_shared_libs = ["liblapack.so", "libblas.so", "libcblas.so", "liblapacke.so"],
    install = True,
    cache_entries = {
        "CMAKE_BUILD_TYPE": "Release",
        "BUILD_TESTING": "ON",
        "BUILD_SHARED_LIBS": "ON",
        "LAPACKE": "ON",
        "CBLAS": "ON",
        "BUILD_DEPRECATED": "ON",
        "CMAKE_INSTALL_PREFIX": "$$EXT_BUILD_DEPS$$",
        "CMAKE_INSTALL_LIBDIR": "lib",
    },
    visibility = ["//visibility:public"],
)

filegroup(
    name = "lapack_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)

filegroup(
    name = "lapack_headers",
    srcs = glob(["**/*.h"]),
    visibility = ["//visibility:public"],
)

cc_library(
    name = "lapack",
    deps = [":lapack_cmake"],
    visibility = ["//visibility:public"],
)
