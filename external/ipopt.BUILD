# ipopt.BUILD
load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")

filegroup(
    name = "ipopt_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
filegroup(
    name = "ipopt_headers",
    srcs = glob(["**/*.h"]),
    visibility = ["//visibility:public"],
)

configure_make(
    name = "ipopt_configure_make",
    lib_source = ":ipopt_sources",
    configure_in_place = False,
    #    autoconf = True,
    #    autogen = True,
    out_include_dir = "include",
    out_lib_dir = "lib",
    configure_options = [
        "--without-hsl",
        "--disable-java",
        "--with-mumps",
        "--with-mumps-cflags=-I$$EXT_BUILD_DEPS$$/include/mumps_seq",
        "--with-mumps-lflags=\"-ldmumps_seq -lmumps_common_seq -lpord_seq -lmpiseq_seq -lesmumps -lscotch -lscotcherr -lmetis -lgfortran\"",
        "--with-asl",
        "--with-asl-cflags=-I$$EXT_BUILD_DEPS$$/include/asl",
        "--with-asl-lflags=-lasl",
        "LDFLAGS=\"-L$$EXT_BUILD_DEPS$$/lib -Wl,-rpath,$$EXT_BUILD_DEPS$$/lib\"",
    ],
    visibility = ["//visibility:public"],
    linkopts = [
        "-L$$EXT_BUILD_DEPS$$/lib",
        "-Wl,-rpath,$$EXT_BUILD_DEPS$$/lib",
        "-lrt",
    ],
    deps = [
        "@com_github_reference_lapack_lapack//:lapack",
        "@com_github_coin_or_tools_thirdparty_mumps//:mumps",
        # Add dependencies here
        #"@mumps//:mumps",
        #"@asl//:asl",
    ],
)

cc_library(
    name = "ipopt",
    hdrs = [":ipopt_configure_make"],
    # We don't need to include the sources here, because this library
    # is a header-only library. CMake is just used for testing etc.
    visibility = ["//visibility:public"],
    includes = ["include"],
)
