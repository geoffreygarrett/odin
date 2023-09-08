# BUILD
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")
load("//:bazel/tools/build_rules.bzl", "create_aliases", "create_config_setting")

# BUILD
config_setting(
    name = "local_brew_gnu_gsl",
    values = {
        "define": "use_local_homebrew_gsl=true",
    },
)

#def gsl_dependency():  # somehow need to check, and raise error for user to make it more explicit
#    # Check if GSL exists at the Homebrew path
#    if paths.exists("/opt/homebrew/opt/gsl"):
#        new_local_repository(
#            name = "org_gnu_gsl",
#            build_file_content = """
#    cc_library(
#    name = "gsl",
#    hdrs = glob(["include/**/*.h"]),
#    includes = ["include"],
#    visibility = ["//visibility:public"],
#    )
#    """,
#            path = "/opt/homebrew/opt/gsl",
#        )

# Creating aliases, from all dependencies defined in the WORKSPACE file.
ALIASES = {
    "eigen": "@com_github_eigen_eigen//:header_lib",
    "autodiff": "@com_github_autodiff_autodiff//:autodiff",
    "glog": "@com_github_google_glog//:glog",
    "cereal": "@com_github_uscilab_cereal//:cereal",
    "gsl": select({
        ":local_brew_gnu_gsl": "@local_brew_gnu_gsl//:gsl",
        "//conditions:default": "@org_gnu_gsl//:gsl",
    }),
    "tbb": "@com_github_oneapi_onetbb//:tbb",
    "benchmark": "@com_google_benchmark//:benchmark",
    "highfive": "@com_github_bluebrain_highfive//:highfive_cmake",
    "zlib": "@net_zlib//:zlib",
    #    "conda_libs": "@conda_prefix//:lib",
    #    "conda_includes": "@conda_prefix//:include",
}

exports_files(
    ["repositories.bzl"],
    visibility = ["//visibility:public"],
)

create_aliases(
    ALIASES,
    visibility = ["//visibility:public"],
)
# TODO: THIS IS USEFUL: TEMPLATING!
#
#load("@bazel_skylib//rules:expand_template.bzl", "expand_template")
#load("@rules_cuda//cuda:defs.bzl", "cuda_library")
#load(":nvbench.bzl", "nvbench_examples")
#
#expand_template(
#    name = "config_cuh",
#    out = "nvbench_generated/nvbench/config.cuh",
#    substitutions = {
#        "#cmakedefine NVBENCH_HAS_NVML": "#define NVBENCH_HAS_NVML 1",
#        "#cmakedefine NVBENCH_HAS_CUPTI": "#define NVBENCH_HAS_CUPTI 1",
#    },
#    template = "nvbench/nvbench/config.cuh.in",
#)

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
CORE_LINKOPTS = select({
    "@platforms//os:osx": [
        #        "-lomp",
    ],
    "//conditions:default": [
        "-fopenmp",  # Add the -fopenmp flag
    ],
})

CORE_COPTS = select({
    "@platforms//os:osx": [
        "-Xpreprocessor",
        "-fopenmp",  # Add the -fopenmp flag
    ],
    "//conditions:default": [
        "-fopenmp",  # Add the -fopenmp flag
        "-fPIC",
        "-std=c++20",
    ],
})

cc_library(
    name = "odin",
    hdrs = [
        ":odin_headers",
    ],
    copts = CORE_COPTS,
    includes = ["include"],
    linkopts = CORE_LINKOPTS,
    visibility = ["//visibility:public"],
    deps = [
        "@com_github_eigen_eigen//:header_lib",
        "@com_github_gabime_spdlog//:spdlog",
        ":gsl",

        #        "@com_github_nvidia_thrust//:thrust",
        "@com_github_oneapi_onetbb//:tbb",
        "@com_github_uscilab_cereal//:cereal",
        #        "@libomp//:omp",
    ],
    #    + select({
    #        "@platforms//os:osx": [
    #            "@libomp//:omp",
    #        ],
    #        "//conditions:default": [
    #        ],
    #    }),
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
