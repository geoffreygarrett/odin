# BUILD
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")
load("//:bazel/tools/build_rules.bzl", "create_aliases", "create_config_setting")

# Creating aliases, from all dependencies defined in the WORKSPACE file.
ALIASES = {
    "eigen": "@com_github_eigen_eigen//:eigen",
    "autodiff": "@com_github_autodiff_autodiff//:autodiff",
    "glog": "@com_github_google_glog//:glog",
    "cereal": "@com_github_uscilab_cereal//:cereal",
    "gsl": "@org_gnu_gsl//:gsl",
    "tbb": "@com_github_oneapi_onetbb//:tbb",
    "benchmark": "@com_google_benchmark//:benchmark",
    "highfive": "@com_github_bluebrain_highfive//:highfive_cmake",
    "zlib": "@net_zlib//:zlib",
    #    "conda_libs": "@conda_prefix//:lib",
    #    "conda_includes": "@conda_prefix//:include",
}

create_aliases(
    ALIASES,
    visibility = ["//visibility:public"],
)

genrule(
    name = "version",
    srcs = ["//:VERSION"],
    outs = ["include/odin/version.hpp"],
    cmd = (
        "VERSION=$$(cat $(location //:VERSION)) && " +
        "major=$$(echo $$VERSION | awk -F. '{print $$1}') && " +
        "minor=$$(echo $$VERSION | awk -F. '{print $$2}') && " +
        "patch=$$(echo $$VERSION | awk -F. '{print $$3}') && " +
        "echo '#define ODIN_VERSION_MAJOR \"'$$major'\"' > $@ && " +
        "echo '#define ODIN_VERSION_MINOR \"'$$minor'\"' >> $@ && " +
        "echo '#define ODIN_VERSION_PATCH \"'$$patch'\"' >> $@ && " +
        "echo '#define ODIN_VERSION \"'$$VERSION'\"' >> $@"
    ),
)

filegroup(
    name = "odin_headers",
    srcs = glob([
        "include/**/*.hpp",
        "include/**/*.h",
        "include/**/*.cuh",
    ]),
)

# Odin Library
cc_library(
    name = "odin",
    hdrs = [
        ":odin_headers",
    ],
    copts = [
        "-std=c++20",
    ],
    includes = ["include"],
    visibility = ["//visibility:public"],
    deps = [
        ":cereal",
        ":eigen",
    ],
)

# Configuration settings
create_config_setting(
    "with_logging",
    "logging=true",
    visibility = ["//visibility:public"],
)

create_config_setting(
    "without_logging",
    "logging=false",
    visibility = ["//visibility:public"],
)

#load("@rules_pkg//pkg:mappings.bzl", "pkg_attributes", "pkg_filegroup", "pkg_files", "pkg_mkdirs", "strip_prefix")
#load("@rules_pkg//pkg:tar.bzl", "pkg_tar")
#load("@rules_pkg//pkg:zip.bzl", "pkg_zip")
#
#pkg_tar(
#    name = "main_tar",
#    srcs = [
#        "README.md",
#        ":odin",
#        ":odin_headers",
#        ":version",
#        "//src:main",
#    ],
#    strip_prefix = "/",
#)
#
#pkg_zip(
#    name = "main_zip",
#    srcs = [
#        "README.md",
#        ":odin",
#        ":odin_headers",
#        ":version",
#        "//src:main",
#    ],
#    strip_prefix = "/",
#)

#load("//:bazel/tools/conda_env.bzl", "conda_environment_setup")
