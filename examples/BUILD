cc_binary(
    name = "dev",
    srcs = ["dev.cpp"],
    deps = ["@//:odin"],
)

cc_binary(
    name = "mcts_tictactoe",
    srcs = ["mcts_tictactoe.cpp"],
    copts = ["-std=c++20"],
    deps = ["@//:odin"],
)

cc_binary(
    name = "example_dev_autodiff",
    srcs = ["example_dev_autodiff.cpp"],
    defines = ["ODIN_AUTODIFF"],
    deps = [
        "@//odin",
        "@//odin:autodiff",
    ],
)

cc_binary(
    name = "example_autodiff_stts",
    srcs = ["example_autodiff_stts.cpp"],
    defines = ["ODIN_AUTODIFF"],
    deps = [
        "@//odin:autodiff",
        "@//odin:eigen",
    ],
)
