# bazel_shtns.BUILD
load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")

filegroup(
    name = "hdf5_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:private"],
    # Private is default, but just to be explicit. Sources would only be
    # part of the cmake target, but not part of the cc_library target.
)

filegroup(
    name = "hdf5_headers",
    srcs = glob(["include/**"]),
    visibility = ["//visibility:public"],
    # We make this public, so we have the option of packaging it
    # with our output.
)

configure_make(
    name = "hdf5_configure_make",
    lib_source = "//:hdf5_sources",
    configure_in_place = True,
    autoconf = False,
    autogen = False,
    linkopts = ["-lz"],
    deps = ["@net_zlib//:zlib"],
    #    copts = ["--enable-openmp"],
    out_include_dir = "include",
    out_lib_dir = "lib",
    out_static_libs = ["libhdf5.a"],
    out_shared_libs = ["libhdf5.so"],
    # Bazel focuses on hermetic, reproducible and deterministic builds. This is why all items such as
    # __DATE__, __TIMESTAMP__, and __TIME__ are replaced with "redacted" to avoid the timestamp
    # changing for each build, the only problem is if these are used In C code, it will not compile.
    # We need to escape the quotes correctly as seen here.
    configure_options = ["CFLAGS=-Dredacted='\\\"redacted\\\"'"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "hdf5",
    deps = [":hdf5_configure_make"],
    visibility = ["//visibility:public"],
    includes = ["include"],
)
