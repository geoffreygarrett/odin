load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

package(
    default_visibility = ["//visibility:public"],
)

#filegroup(
#    name = "eigen_hdrs",
#    srcs = glob(["Eigen/**/*.h"]),
#    visibility = ["//visibility:public"],
#)
#
#cc_import(
#    name = "eigen",
#    hdrs = ":eigen_hdrs",
#    visibility = ["//visibility:public"],
#)

#cc_library(
#    name = "eigen",
#    hdrs = glob(["Eigen/**/*.h"]),
#    includes = ["include"],
#    copts = ["-I./"]
#)

cc_library(
    name = 'eigen',
    srcs = [],
    includes = ['.'],
    hdrs = glob(['Eigen/**']),
    visibility = ['//visibility:public'],
)

#cmake(
#    name = "build_eigen",
#    lib_source = "//:eigen_sources",
#    out_include_dir = "include",
#    out_headers_only = True, # Flag variable to indicate that the library produces only headers
#    install = True,
#    cache_entries = {
##        "CMAKE_INSTALL_DATADIR": "lib",
##        "INCLUDE_INSTALL_DIR": "include/eigen3",
#    }
#)
#
#filegroup(
#     name = "eigen_sources",
#     srcs = glob(["**"]),
#     visibility = ["//visibility:public"],
# )