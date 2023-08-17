load("@bazel_tools//tools/cpp:toolchain_utils.bzl", "find_cpp_toolchain")

cc_library(
    name = "xsimd",
    hdrs = glob([
        "include/xsimd/arch/*.hpp",
        "include/xsimd/config/*.hpp",
        "include/xsimd/memory/*.hpp",
        "include/xsimd/types/*.hpp",
        "include/xsimd/*.hpp",
        "include/**/*.hpp",
    ]),
    includes = ["include"],
    visibility = ["//visibility:public"],
)

config_setting(
    name = "enable_xtl_complex",
    define_values = {"ENABLE_XTL_COMPLEX": "1"},
)

alias(
    name = "xsimd_config",
    actual = select({
        ":enable_xtl_complex": ":xsimd_with_xtl_complex",
        "//conditions:default": ":xsimd",
    }),
)

cc_library(
    name = "xsimd_with_xtl_complex",
    deps = [
        ":xsimd",
        "@xtl//:xtl",
    ],
    defines = ["XSIMD_ENABLE_XTL_COMPLEX=1"],
)

cc_test(
    name = "xsimd_tests",
    srcs = glob(["test/**/*.cpp"]),
    deps = [":xsimd_config"],
)

cc_binary(
    name = "xsimd_benchmarks",
    srcs = glob(["benchmark/**/*.cpp"]),
    deps = [":xsimd_config"],
)

cc_binary(
    name = "xsimd_examples",
    srcs = glob(["examples/**/*.cpp"]),
    deps = [":xsimd_config"],
)
