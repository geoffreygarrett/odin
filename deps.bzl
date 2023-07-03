load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:utils.bzl", "maybe")

def odin_dependencies():
    #    maybe(
    #        http_archive,
    #        name = "rules_pkg",
    #        sha256 = "360c23a88ceaf7f051abc99e2e6048cf7fe5d9af792690576554a88b2013612d",
    #        strip_prefix = "rules_pkg-0.9.1",
    #        urls = ["https://github.com/bazelbuild/rules_pkg/archive/0.9.1.tar.gz"],
    #    )
    #    maybe(
    #        http_archive,
    #        name = "rules_foreign_cc",
    #        sha256 = "6041f1374ff32ba711564374ad8e007aef77f71561a7ce784123b9b4b88614fc",
    #        strip_prefix = "rules_foreign_cc-0.8.0",
    #        urls = ["https://github.com/bazelbuild/rules_foreign_cc/archive/0.8.0.tar.gz"],
    #    )

    maybe(
        http_archive,
        name = "rules_pkg",
        sha256 = "360c23a88ceaf7f051abc99e2e6048cf7fe5d9af792690576554a88b2013612d",
        strip_prefix = "rules_pkg-0.9.1",
        urls = ["https://github.com/bazelbuild/rules_pkg/archive/0.9.1.tar.gz"],
    )
    maybe(
        http_archive,
        name = "rules_foreign_cc",
        sha256 = "6041f1374ff32ba711564374ad8e007aef77f71561a7ce784123b9b4b88614fc",
        strip_prefix = "rules_foreign_cc-0.8.0",
        urls = ["https://github.com/bazelbuild/rules_foreign_cc/archive/0.8.0.tar.gz"],
    )
    maybe(
        http_archive,
        name = "com_github_gflags_gflags",
        sha256 = "34af2f15cf7367513b352bdcd2493ab14ce43692d2dcd9dfc499492966c64dcf",
        strip_prefix = "gflags-2.2.2",
        urls = ["https://github.com/gflags/gflags/archive/v2.2.2.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_google_glog",
        sha256 = "122fb6b712808ef43fbf80f75c52a21c9760683dae470154f02bddfc61135022",
        strip_prefix = "glog-0.6.0",
        urls = ["https://github.com/google/glog/archive/v0.6.0.zip"],
    )

    maybe(
        http_archive,
        name = "com_github_eigen_eigen",
        build_file = "//:external/eigen.BUILD",
        sha256 = "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72",
        strip_prefix = "eigen-3.4.0",
        urls = ["https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_autodiff_autodiff",
        build_file = "//:external/bazel_autodiff.BUILD",
        sha256 = "ce87b642272a84f3dcce5463d0f97ef30b15e5a5645898061c06c16a8d61bf22",
        strip_prefix = "autodiff-1.0.3",
        urls = ["https://github.com/autodiff/autodiff/archive/v1.0.3.zip"],
    )

    maybe(
        http_archive,
        name = "com_github_uscilab_cereal",
        build_file = "//:external/cereal.BUILD",
        sha256 = "e72c3fa8fe3d531247773e346e6824a4744cc6472a25cf9b30599cd52146e2ae",
        strip_prefix = "cereal-1.3.2",
        urls = ["https://github.com/USCiLab/cereal/archive/v1.3.2.zip"],
    )

    maybe(
        http_archive,
        name = "com_google_googletest",
        sha256 = "ffa17fbc5953900994e2deec164bb8949879ea09b411e07f215bfbb1f87f4632",
        strip_prefix = "googletest-1.13.0",
        urls = ["https://github.com/google/googletest/archive/v1.13.0.zip"],
    )

    maybe(
        http_archive,
        name = "com_google_benchmark",
        sha256 = "927475805a19b24f9c67dd9765bb4dcc8b10fd7f0e616905cc4fee406bed81a7",
        strip_prefix = "benchmark-1.8.0",
        urls = ["https://github.com/google/benchmark/archive/v1.8.0.zip"],
    )

    maybe(
        http_archive,
        name = "org_gnu_gsl",
        build_file = "//:external/gsl.BUILD",
        sha256 = "dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b",
        strip_prefix = "gsl-2.7.1",
        urls = ["https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_oneapi_onetbb",
        sha256 = "fcebb93cb9f7e882f62cd351b1c093dbefdcae04b616227dc716b0a5efa9e8ab",
        strip_prefix = "oneTBB-2021.9.0",
        urls = ["https://github.com/oneapi-src/oneTBB/archive/v2021.9.0.zip"],
    )

    maybe(
        http_archive,
        name = "com_github_bluebrain_highfive",
        build_file = "//:external/highfive.BUILD",
        sha256 = "029fb6b13fe6ef8098e879e95a55f501f34760d2f9aaa2bdfd513bba1058276b",
        strip_prefix = "HighFive-2.7.1",
        urls = ["https://github.com/BlueBrain/HighFive/archive/v2.7.1.zip"],
    )

    maybe(
        http_archive,
        name = "net_zlib",
        build_file = "//:external/zlib.BUILD",
        sha256 = "b3a24de97a8fdbc835b9833169501030b8977031bcb54b3b3ac13740f846ab30",
        strip_prefix = "zlib-1.2.13",
        urls = ["https://zlib.net/zlib-1.2.13.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_hdfgroup_hdf5",
        build_file = "//:external/hdf5.BUILD",
        sha256 = "5fbebeb7162d336ed06f119ba32f336f2a67c0a492c7802e48f16d4cfb8b0fbf",
        strip_prefix = "hdf5-hdf5-1_14_1-2",
        urls = ["https://github.com/HDFGroup/hdf5/archive/hdf5-1_14_1-2.zip"],
    )
