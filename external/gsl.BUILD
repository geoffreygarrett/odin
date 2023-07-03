load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")

configure_make(
    name = "gsl",
    lib_source = "//:gsl_sources",
    out_include_dir = "include",
    out_lib_dir = "lib",
    out_static_libs = ["libgsl.a", "libgslcblas.a"],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "gsl_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
