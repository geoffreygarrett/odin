#load("@rules_foreign_cc//tools/build_defs:cmake.bzl", "cmake")
#load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

config_setting(
    name = "darwin",
    constraint_values = ["@platforms//os:mac"],
)

config_setting(
    name = "cereal_master_project",
    values = {"define": "cereal_master_project=1"},
)

filegroup(
    name = "cereal_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:private"],
    # Private is default, but just to be explicit. Sources would only be
    # part of the cmake target, but not part of the cc_library target.
)

filegroup(
    name = "cereal_headers",
    srcs = glob(["include/cereal/**"]),
    visibility = ["//visibility:public"],
    # We make this public, so we have the option of packaging it
    # with our output.
)

#genrule(
#    name = "generate_cmake_cache_entries",
#    outs = ["CMakeCacheEntries.txt"],
#    cmd = "\n".join([
#        "echo {}={} >> $@".format(key, SELECT_DICT[key][build_setting])
#        for key, build_setting in SELECT_DICT.items()
#    ]),
#)

# This will be used later to run the respective libraries tests in our
# build configuration and environments.
#cmake(
#    name = "cmake",
#    lib_source = ":cereal_sources",
#    out_include_dir = "include",
#    cache_entries = {
#        "CMAKE_BUILD_TYPE": "Release",
#        "SKIP_PORTABILITY_TEST": str(select({
#            "//conditions:default": "ON",
#            ":darwin": "ON",
#        })),
#        "BUILD_DOC": str(select({
#            "//conditions:default": "OFF",
#            ":cereal_master_project": "ON",
#        })),
#        "BUILD_SANDBOX": str(select({
#            "//conditions:default": "OFF",
#            ":cereal_master_project": "ON",
#        })),
#        "SKIP_PERFORMANCE_COMPARISON": "OFF",
#        "JUST_INSTALL_CEREAL": "OFF",
#        "THREAD_SAFE": "OFF",
#        "WITH_WERROR": "ON",
#        "CLANG_USE_LIBCPP": "OFF",
#        "CMAKE_CXX_STANDARD": "11",
#        "CEREAL_INSTALL": str(select({
#            "//conditions:default": "OFF",
#            ":cereal_master_project": "ON",
#        })),
#        "BUILD_TESTS": str(select({
#            "//conditions:default": "OFF",
#            ":cereal_master_project": "ON",
#        })),
#    },
#    visibility = ["//visibility:public"],
#    deps = [
#        "@com_google_googletest//:gtest_main",
#        "@com_google_googletest//:gtest",
#    ],
#)

cc_library(
    name = "cereal",
    hdrs = [":cereal_headers"],
    # We don't need to include the sources here, because this library
    # is a header-only library. CMake is just used for testing etc.
    visibility = ["//visibility:public"],
    includes = ["include"],
)
