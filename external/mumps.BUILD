# BUILD file
load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")

configure_make(
    name = "mumps_configure_make",
    lib_source = "//:mumps_sources",  # Ensure this points to the correct SHTns library source files
    configure_in_place = True,
    configure_options = [
        "--without-hsl",
        "--disable-java",
        "--with-mumps",
        "--with-mumps-cflags=-I$$EXT_BUILD_DEPS$$/include/mumps_seq",
        "--with-mumps-lflags=\"-ldmumps_seq -lmumps_common_seq -lpord_seq -lmpiseq_seq -lesmumps -lscotch -lscotcherr -lmetis -lgfortran\"",
        "LDFLAGS=\"-L$$EXT_BUILD_DEPS$$/lib -Wl,-rpath,$$EXT_BUILD_DEPS$$/lib\"",
    ],
    out_include_dir = "include",
    out_lib_dir = "lib",
    visibility = ["//visibility:public"],
)

filegroup(
    name = "mumps_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)

cc_library(
    name = "mumps",
    srcs = ["//:mumps_configure_make"],
    #    hdrs = glob(["include/**"]),
    #    includes = ["include"],
    visibility = ["//visibility:public"],
)
