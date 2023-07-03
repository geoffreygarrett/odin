#load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")
## https://bazelbuild.github.io/rules_foreign_cc/0.5.1/cmake.html
#package(
#    default_visibility = ["//visibility:public"],
#)
#
#
#genrule(
#    name = "rename_autodiff",
#    srcs = glob(["**"]),
#    outs = ["autodiff_renamed"],
#    cmd = "mv $(location autodiff) $(location autodiff_renamed)",
#    visibility = ["//visibility:public"],
#)
#
#cmake(
#    name = "autodiff",
#    lib_source = "//:autodiff_sources", # Label with source code to build. Typically a `filegroup` for the source of remote repository. Mandatory.
#    out_include_dir = "autodiff", # Optional name of the output subdirectory with the header files, defaults to 'include'.
#    out_headers_only = True, # Flag variable to indicate that the library produces only headers
#    install = True, # If True, the cmake --install command will be performed after a build
#    cache_entries = {
#        "AUTODIFF_BUILD_TESTS": "OFF",
#        "AUTODIFF_BUILD_PYTHON": "OFF",
#        "AUTODIFF_BUILD_EXAMPLES": "OFF",
#        "AUTODIFF_BUILD_DOCS": "OFF",
#        "Eigen3_DIR": "$$EXT_BUILD_DEPS$$/eigen/share/eigen3/cmake",
#        },
#    deps = [
#        "@eigen//:eigen",
##        ":rename_autodiff",  # add this line
#    ],
#)

# Odin Library
cc_library(
    name = "autodiff",
    hdrs = glob(["autodiff/**/*.hpp"]),
#    includes = ["include"],
    visibility = ["//visibility:public"],
)

#filegroup(
#    name = "autodiff_sources",
#    srcs = glob(["**"]),
#    visibility = ["//visibility:public"],
#)
