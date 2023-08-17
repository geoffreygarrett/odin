load("@bazel_skylib//rules:common_settings.bzl", "string_flag")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")
#load("@//:external/_thrust/defs.bzl", "get_cmake_cache_entries", "get_cmake_defines")

### Flags
#string_flag(
#    name = "thrust_enable_header_testing",
#    build_setting_default = "ON",
#)
#
#string_flag(
#    name = "thrust_enable_testing",
#    build_setting_default = "ON",
#)
#
#string_flag(
#    name = "thrust_enable_examples",
#    build_setting_default = "ON",
#)
#
#string_flag(
#    name = "thrust_enable_benchmarks",
#    build_setting_default = "OFF",
#)
#
#string_flag(
#    name = "thrust_include_cub_cmake",
#    build_setting_default = "OFF",
#)

## File-groups
filegroup(
    name = "thrust_all_headers",
    srcs = glob(["**/*.hpp", "**/*.h"]),
    visibility = ["//visibility:public"],
)

filegroup(
    name = "thrust_all_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:private"],
)

# Defines
#DEFINES = [
#    "THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP",
#    "THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB",
#]

# Cache Entries
CACHE_ENTRIES = {
    "THRUST_HOST_SYSTEM": "CPP",
    "THRUST_DEVICE_SYSTEM": "CPP",
    "CMAKE_CXX_STANDARD": "14",
    "THRUST_ENABLE_HEADER_TESTING": "OFF",
    "THRUST_ENABLE_TESTING": "OFF",
    "THRUST_ENABLE_EXAMPLES": "OFF",
    "THRUST_ENABLE_MULTICONFIG": "OFF",
    "THRUST_ENABLE_EXAMPLE_FILECHECK": "OFF",
    "THRUST_ENABLE_INSTALL_RULES": "ON",
    "THRUST_INCLUDE_CUB_CMAKE": "OFF",
    "THRUST_INSTALL_CUB_HEADERS": "OFF",
}

cc_library(
    name = "libcudacxx",
    hdrs = glob(["dependencies/libcudacxx/*.h"]),
    includes = ["dependencies/libcudacxx/"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "cub",
    hdrs = glob(["dependencies/cub/*.h"]),
    includes = ["dependencies/cub/"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "thrust",
    hdrs = glob(["**/*.h"], exclude = ["dependencies/**"]),
    includes = ["."],
    visibility = ["//visibility:public"],
    defines = [
        "THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPU",
        "THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPU",
    ],
)

#cmake(
#    name = "thrust",
#    lib_source = ":thrust_all_sources",
#    out_headers_only = True,
#    install = True,
#    defines = [
#        "THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPU",
#        "THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB",
#    ],
#    #    defines = DEFINES,
#    #    out_include_dir = "include",
#    cache_entries = CACHE_ENTRIES,
#    #    cache_entries = {
#    #        #GENERIC_CACHE_ENTRIES = {
#    #        "CMAKE_CXX_STANDARD": "17",
#    #        "THRUST_ENABLE_HEADER_TESTING": "OFF",
#    #        "THRUST_ENABLE_TESTING": "OFF",
#    #        "THRUST_ENABLE_EXAMPLES": "OFF",
#    #        "THRUST_ENABLE_MULTICONFIG": "OFF",
#    #        "THRUST_ENABLE_EXAMPLE_FILECHECK": "OFF",
#    #        "THRUST_ENABLE_INSTALL_RULES": "OFF",
#    #        #        }
#    #
#    #        #SINGLE_CONFIG_CMAKE_OPTIONS = {
#    #        "THRUST_HOST_SYSTEM": "CPP",
#    #        "THRUST_DEVICE_SYSTEM": "TBB",
#    #        #        "THRUST_CPP_DIALECT": "17",
#    #        #}
#    #
#    #        #CUDA_SPECIFIC_CMAKE_OPTIONS = {
#    #        #        "THRUST_INCLUDE_CUB_CMAKE": "OFF",
#    #        #        "THRUST_INSTALL_CUB_HEADERS": "OFF",
#    #        #        "THRUST_ENABLE_COMPUTE_XX": "OFF",
#    #        #        "THRUST_ENABLE_COMPUTE_FUTURE": "OFF",
#    #        #        "THRUST_DISABLE_ARCH_BY_DEFAULT": "OFF",
#    #        #        "THRUST_ENABLE_TESTS_WITH_RDC": "OFF",
#    #        #        "THRUST_ENABLE_EXAMPLES_WITH_RDC": "OFF",
#    #        #}
#    #    },
#    visibility = ["//visibility:public"],
#    deps = [
#        "@com_github_oneapi_onetbb//:tbb",
#    ],
#)

#alias(
#    name = "thrust",
#    actual = ":thrust_cmake",
#)

#cc_library(
#    name = "thrust",
#    deps = [":thrust_cmake"],
#    visibility = ["//visibility:public"],
#    includes = ["."],
#)
