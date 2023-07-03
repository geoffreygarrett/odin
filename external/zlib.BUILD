# Bazel build file for zlib compression library

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "zlib",
    hdrs = glob(["**/*.h"]),
    includes = ["."],
    srcs = [
        "adler32.c",
        "compress.c",
        "crc32.c",
        "crc32.h",
        "deflate.c",
        "deflate.h",
        "gzclose.c",
        "gzguts.h",
        "gzlib.c",
        "gzread.c",
        "gzwrite.c",
        "infback.c",
        "inffast.c",
        "inffast.h",
        "inffixed.h",
        "inflate.c",
        "inflate.h",
        "inftrees.c",
        "inftrees.h",
        "trees.c",
        "trees.h",
        "uncompr.c",
        "zconf.h",
        "zlib.h",
        "zutil.c",
        "zutil.h",
    ],
    copts = [
        "-Wall",
        "-Wextra",
        "-Wno-sign-compare",
        "-Wno-unused-parameter",
        "-Wno-implicit-function-declaration",
    ],
)

cc_library(
    name = "minizip",
    srcs = [
        "contrib/minizip/crypt.h",
        "contrib/minizip/ioapi.c",
        "contrib/minizip/ioapi.h",
        "contrib/minizip/unzip.c",
        "contrib/minizip/unzip.h",
        "contrib/minizip/zip.c",
        "contrib/minizip/zip.h",
    ],
    deps = [
        ":zlib",
    ],
    copts = [
        "-Wall",
        "-Wextra",
        "-Wno-sign-compare",
        "-Wno-unused-parameter",
        "-Wno-implicit-function-declaration",
        "-Wno-parentheses-equality",
    ],
)
