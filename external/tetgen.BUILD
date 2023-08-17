load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

# Rule to replace the string in the source file
genrule(
    name = "replace_tetgen",
    srcs = ["tetgen.cxx"],
    outs = ["tetgen_mod.cxx"],
    cmd = 'cat $(SRCS) | sed "s/#include \\"tetgen.h\\"/#include \\"tetgen.h\\"\\n#define printf(fmt, ...) (0)\\n/" > $(OUTS)',
)

cc_library(
    name = "tetgen",
    srcs = [
        "predicates.cxx",
        "tetgen_mod.cxx",
    ],
    hdrs = ["tetgen.h"],
    copts = [
        "-fPIC",
        "-Wunused-variable",
        "-Wsign-compare",
        "-w",  # to suppress warnings
    ],
    defines = ["TETLIBRARY"],
    includes = ["@tetgen//"],
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "tetgen_exec",
    deps = [
        ":tetgen",
    ],
)
