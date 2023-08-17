load("@bazel_skylib//rules:common_settings.bzl", "bool_flag", "string_flag")

string_flag(
    name = "polyhedral_gravity_parallelization",
    build_setting_default = "CPP",
)
string_flag(
    name = "logging_level",
    build_setting_default = "2",
)

#bool_flag(
#    name = "use_local_tbb",
#    build_setting_default = True,
#)
bool_flag(
    name = "build_polyhedral_gravity_docs",
    build_setting_default = False,
)

bool_flag(
    name = "build_polyhedral_gravity_executable",
    build_setting_default = True,
)

bool_flag(
    name = "build_polyhedral_gravity_library",
    build_setting_default = True,
)

bool_flag(
    name = "build_polyhedral_gravity_python_interface",
    build_setting_default = True,
)

bool_flag(
    name = "build_polyhedral_gravity_tests",
    build_setting_default = True,
)

cc_library(
    name = "polyhedral_gravity",
    srcs = glob(["src/**/*.cpp"]),
    hdrs = glob(["src/**/*.h"]),
    includes = ["./src"],
    copts = [
        "-std=c++17",
        "-DPOLYHEDRAL_GRAVITY_PARALLELIZATION=CPP",
        "-DLOGGING_LEVEL=2",
    ],
    deps = [
        "@com_github_jbeder_yaml_cpp//:yaml-cpp",
        "@com_github_nvidia_thrust//:thrust_cpp_tbb",
        "@com_github_gabime_spdlog//:spdlog",
        "@com_github_xtensor_stack_xsimd//:xsimd",
        "@com_github_libigl_tetgen//:tetgen",
        "@com_github_oneapi_onetbb//:tbb",
        "@com_google_googletest//:gtest",
    ],
    defines = [
        #        set(POLYHEDRAL_GRAVITY_PARALLELIZATION "CPP" CACHE STRING "Host parallelization chosen by the user
        #         (CPP= Serial, OMP = OpenMP, TBB = Intel Threading Building Blocks")
        "POLYHEDRAL_GRAVITY_PARALLELIZATION=CPP",
        #        Enforce to use an already installed tbb library instead of compiling from source
        "USE_LOCAL_TBB=ON",
        #        set(LOGGING_LEVEL "2" CACHE STRING "Set the Logging level, default (INFO=2), available options:
        #        TRACE=0, DEBUG=1, INFO=2, WARN=3, ERROR=4, CRITICAL=5, OFF=6")
        "LOGGING_LEVEL=2",
        "BUILD_POLYHEDRAL_GRAVITY_DOCS=OFF",
        "BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE=OFF",
        "BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE=OFF",
        "BUILD_POLYHEDRAL_GRAVITY_TESTS=OFF",
        "BUILD_POLYHEDRAL_GRAVITY_LIBRARY=OFF",
    ],

    #    defines = [
    #        "USE_LOCAL_TBB=" + select({
    #            ":use_local_tbb": "ON",
    #            "//conditions:default": "OFF",
    #        }),
    #        "BUILD_POLYHEDRAL_GRAVITY_DOCS=" + select({
    #            ":build_polyhedral_gravity_docs": "ON",
    #            "//conditions:default": "OFF",
    #        }),
    #        "BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE=" + select({
    #            ":build_polyhedral_gravity_executable": "ON",
    #            "//conditions:default": "ON",
    #        }),
    #        "BUILD_POLYHEDRAL_GRAVITY_LIBRARY=" + select({
    #            ":build_polyhedral_gravity_library": "ON",
    #            "//conditions:default": "ON",
    #        }),
    #        "BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE=" + select({
    #            ":build_polyhedral_gravity_python_interface": "ON",
    #            "//conditions:default": "ON",
    #        }),
    #        "BUILD_POLYHEDRAL_GRAVITY_TESTS=" + select({
    #            ":build_polyhedral_gravity_tests": "ON",
    #            "//conditions:default": "ON",
    #        }),
    #    ],
    visibility = ["//visibility:public"],
)
