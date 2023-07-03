#
## bazel_fftw.BUILD
#load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")
#
#configure_make(
#    name = "fftw",
#    lib_source = "//:fftw_sources",
#    out_include_dir = "include",
#    out_lib_dir = "lib",
#    out_static_libs = ["libfftw3.a"],  # Replace with the actual static library name if different
#    #    out_shared_libs = ["libfftw3.so"],  # Replace with the actual shared library name if different
#    visibility = ["//visibility:public"],
#    configure_options = [
#        #        "--enable-shared",
#        #        "--enable-long-double",
#        #        "--enable-quad-precision",
#        #        "--enable-sse2",
#        #        "--enable-avx",
#        #        "--enable-avx2",
#        #        "--enable-avx512",
#        #        "--enable-fma",
#        "--enable-openmp",
#        #        "--enable-threads",
#    ],
#)
#
#filegroup(
#    name = "fftw_sources",
#    srcs = glob(["**"]),
#    visibility = ["//visibility:public"],
#)

# FFT library
#http_archive(
#    name = "org_fftw_fftw",
#    build_file = "bazel_fftw.BUILD",
#    sha256 = "56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467",
#    strip_prefix = "fftw-3.3.10",
#    urls = ["https://www.fftw.org/fftw-3.3.10.tar.gz"],
#)

## Spherical Harmonics Transform library
#http_archive(
#    name = "org_bitbucket_shtns",
#    build_file = "bazel_shtns.BUILD",
#    sha256 = "9f1f46e18c6b346f35c6764f0503678412412719dccc964501e1480dd3087861",
#    strip_prefix = "shtns",
#    urls = ["https://bitbucket.org/nschaeff/shtns/downloads/shtns-3.6.tar.gz"],
#)

# bazel_shtns.BUILD
load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")

filegroup(
    name = "shtns_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)

configure_make(
    name = "shtns",
    lib_source = "//:shtns_sources",
    configure_in_place = True,
    autoconf = False,
    autogen = False,
    copts = ["--enable-openmp"],
    out_include_dir = "include",
    out_lib_dir = "lib",
    out_static_libs = ["libshtns.a"],
    # Bazel focuses on hermetic and to have builds reproducible and deterministic. This is why all items such as
    # __DATE__, __TIMESTAMP__, and __TIME__ are replaced with "redacted" to avoid the timestamp
    # changing for each build, the only problem is if these are used In C code, it will not compile.
    # We need to escape the quotes correctly as seen here.
    configure_options = ["CFLAGS=-Dredacted='\\\"redacted\\\"'"],
    visibility = ["//visibility:public"],
    deps = ["@org_fftw_fftw//:fftw"],  # Add the FFTW library as a dependency
    linkopts = ["-L@org_fftw_fftw//:fftw"],
)

#configure_make(
#    name = "shtns",
#    lib_source = "//:shtns_sources",
#    configure_in_place = True,
#    configure_options = [
#        #        "--enable-openmp",
#        #        "--enable-shared",
#        #        "--enable-avx",
#    ],
#    out_include_dir = "include",
#    out_lib_dir = "lib",
#    out_static_libs = ["libshtns.a"],  # Replace with the actual static library name if different
#    visibility = ["//visibility:public"],
#    deps = ["@org_fftw_fftw//:fftw"],  # Add the FFTW library as a dependency
#    linkopts = ["-L@org_fftw_fftw//:fftw"],
#    #    install_prefix = "$$INSTALLDIR$$",  # Replace with the actual install prefix if different
#    targets = [
#        "install",
#    ],
#)
