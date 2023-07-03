load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")

filegroup(
    name = "hdf5_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:private"],
)

filegroup(
    name = "hdf5_headers",
    srcs = glob(["include/**"]),
    visibility = ["//visibility:public"],
)

configure_make(
    name = "hdf5_configure_make",
    lib_source = "//:hdf5_sources",
    configure_in_place = True,
    autoconf = False,
    autogen = False,
    deps = ["@net_zlib//:zlib"],
    out_include_dir = "include",
    out_lib_dir = "lib",
    out_static_libs = ["libhdf5.a"],
    out_shared_libs = ["libhdf5.so"],
    linkopts = ["-lz"],
    visibility = ["//visibility:public"],
    configure_options = [
        # Bazel focuses on hermetic, reproducible and deterministic builds. This is why all items such as
        # __DATE__, __TIMESTAMP__, and __TIME__ are replaced with "redacted" to avoid the timestamp
        # changing for each build, the only problem is if these are used In C code, it will not compile.
        # We need to escape the quotes correctly as seen here.
        "CFLAGS=-Dredacted='\\\"redacted\\\"'",
        "--with-pic",
        "--with-zlib=${PREFIX}",
        "--with-szlib=${PREFIX}",
        "--with-pthread=yes",
        "--enable-cxx",
        "--enable-fortran",
        "--with-default-plugindir=${PREFIX}/lib/hdf5/plugin",
        "--enable-threadsafe",
        "--enable-build-mode=production",
        "--enable-unsupported",
        "--enable-hl-tools=yes",
        "--enable-using-memchecker",
        "--enable-static=no",
        "--enable-ros3-vfd",
        "--enable-direct-vfd",
        "--enable-parallel",
        "--enable-tests=no",
    ],
    copts = select({
        "@bazel_tools//src/conditions:linux_x86_64": [
            "-Wl,--no-as-needed -Wl,--disable-new-dtags",
        ],
        "//conditions:default": [],
    }),
    build_script = select({
        "@bazel_tools//src/conditions:linux_x86_64": ["./configure --prefix=${PREFIX} --host=${HOST} --build=${BUILD}"],
        "//conditions:default": ["./configure --prefix=${PREFIX}"],
    }),
)

cc_library(
    name = "hdf5",
    deps = [":hdf5_configure_make"],
    visibility = ["//visibility:public"],
    includes = ["include"],
)
