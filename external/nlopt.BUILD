load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

cmake(
    name = "nlopt_cmake",
    cache_entries = {
        #        "CMAKE_PREFIX_PATH": "$$EXT_BUILD_DEPS$$",
        #        "CMAKE_INSTALL_LIBDIR": "lib",
        "CMAKE_BUILD_TYPE": "Release",
        "NLOPT_GUILE": "OFF",
        "NLOPT_MATLAB": "OFF",
        "NLOPT_OCTAVE": "OFF",
    },
    lib_source = ":nlopt_sources",
    out_shared_libs = ["libnlopt.so"],
    out_include_dir = "include",
    install = True,
    visibility = ["//visibility:public"],
)

filegroup(
    name = "nlopt_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)

filegroup(
    name = "nlopt_headers",
    srcs = glob(["**/*./h"]),
    visibility = ["//visibility:public"],
)

cc_library(
    name = "nlopt",
    deps = [":nlopt_cmake"],
    visibility = ["//visibility:public"],
)
