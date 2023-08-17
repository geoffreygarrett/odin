load("@bazel_skylib//rules:common_settings.bzl", "bool_flag", "string_flag")

bool_flag(
    name = "enable_thrust_defines",
    build_setting_default = True,
    visibility = ["//visibility:public"],
)

# Config setting that will be true when enable_thrust_defines is False
config_setting(
    name = "disable_defines",
    flag_values = {":enable_thrust_defines": "False"},
    visibility = ["//visibility:public"],
)

string_flag(
    name = "thrust_host_system",
    build_setting_default = "CPP",
    values = [
        "CPP",
        "OMP",
        "TBB",
    ],
    visibility = ["//visibility:public"],
)

string_flag(
    name = "thrust_device_system",
    build_setting_default = "TBB",
    values = [
        "CPP",
        "OMP",
        "TBB",
        "CUDA",
    ],
    visibility = ["//visibility:public"],
)

config_setting(
    name = "use_cpp_host",
    flag_values = {":thrust_host_system": "CPP"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "use_omp_host",
    flag_values = {":thrust_host_system": "OMP"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "use_tbb_host",
    flag_values = {":thrust_host_system": "TBB"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "use_cpp_device",
    flag_values = {":thrust_device_system": "CPP"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "use_omp_device",
    flag_values = {":thrust_device_system": "OMP"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "use_tbb_device",
    flag_values = {":thrust_device_system": "TBB"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "use_cuda_device",
    flag_values = {":thrust_device_system": "CUDA"},
    visibility = ["//visibility:public"],
)

cc_library(
    name = "libcudacxx",
    hdrs = glob([
        "dependencies/libcudacxx/include/**/*",
    ]),
    includes = ["dependencies/libcudacxx/include"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "cub",
    hdrs = glob([
        "dependencies/cub/**/*.h",
        "dependencies/cub/**/*.hpp",
        "dependencies/cub/**/*.cuh",
    ]),
    includes = ["dependencies/cub/"],
    visibility = ["//visibility:public"],
)

config_setting(
    name = "tbb_host_and_device",
    flag_values = {
        ":thrust_device_system": "TBB",
        ":thrust_host_system": "TBB",
    },
    visibility = ["//visibility:public"],
)

config_setting(
    name = "tbb_host_only",
    flag_values = {
        ":thrust_host_system": "TBB",
    },
    visibility = ["//visibility:public"],
)

config_setting(
    name = "tbb_device_only",
    flag_values = {
        ":thrust_device_system": "TBB",
    },
    visibility = ["//visibility:public"],
)

cc_library(
    name = "thrust",
    hdrs = glob(["**/*.h", "**/*.inl"]),
    includes = ["."],
    visibility = ["//visibility:public"],
    deps = [
        ":libcudacxx",
    ] + select({
        ":tbb_host_and_device": ["@com_github_oneapi_onetbb//:tbb"],
        ":tbb_host_only": ["@com_github_oneapi_onetbb//:tbb"],
        ":tbb_device_only": ["@com_github_oneapi_onetbb//:tbb"],
        "//conditions:default": [],
    }) + select({
        ":use_cuda_device": [":cub", "@local_cuda//:cuda_headers"],  #  "@local_cuda//:cudart",
        "//conditions:default": [],
    }),
    defines = select({
        ":disable_defines": [],
        ":use_cpp_host": ["THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP"],
        ":use_omp_host": ["THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP"],
        ":use_tbb_host": ["THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_TBB"],
        "//conditions:default": ["THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP"],
    }) + select({
        ":disable_defines": [],
        ":use_cpp_device": ["THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP"],
        ":use_omp_device": ["THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP"],
        ":use_tbb_device": ["THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB"],
        ":use_cuda_device": ["THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA"],
        "//conditions:default": ["THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB"],
    }),
)

cc_library(
    name = "thrust_cpp_tbb",
    hdrs = glob(["**/*.h", "**/*.inl"]),
    includes = ["."],
    visibility = ["//visibility:public"],
    deps = [
        ":libcudacxx",
        "@com_github_oneapi_onetbb//:tbb",
    ],
    defines = [
        "THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_TBB",
        "THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB",
    ],
)
