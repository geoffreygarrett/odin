load("@rules_metal//:metal.bzl", "metal_binary", "metal_library")

cc_library(
    name = "metal_adder",
    srcs = ["metal_adder.cpp"],
    hdrs = ["metal_adder.hpp"],
    data = [
        ":default",
    ],
    includes = ["."],
    deps = [
        "@developer_apple_com_metal_cpp//:metal_cpp",
    ],
)

# Produces the final executable
cc_binary(
    name = "app",
    srcs = [
        "main1.cpp",
    ],
    deps = [
        ":metal_adder",
    ],
)

# Produces shaders.metallib
metal_binary(
    name = "default",
    srcs = [
        "add.metal",
    ],
)
