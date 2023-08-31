load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:utils.bzl", "maybe")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_file")
load("@bor//odin:bazel/rules/gsl.bzl", "gsl_repository")

#==============================================================================
#                              Bazel Dependencies
#    Description: Dependencies required for the Bazel build system
#==============================================================================
BAZEL_DEPS = [
    {
        "type": http_archive,
        "name": "bazel_skylib",
        "urls": [
            "https://mirror.bazel.build/github.com/bazelbuild/bazel-skylib/releases/download/1.4.2/bazel-skylib-1.4.2.tar.gz",
            "https://github.com/bazelbuild/bazel-skylib/releases/download/1.4.2/bazel-skylib-1.4.2.tar.gz",
        ],
        "sha256": "66ffd9315665bfaafc96b52278f57c7e2dd09f5ede279ea6d39b2be471e7e3aa",
    },
]

def bazel_dependencies():
    for dep in BAZEL_DEPS:
        maybe(
            dep["type"],
            name = dep["name"],
            urls = dep["urls"],
            sha256 = dep["sha256"],
        )

#==============================================================================
#                         Polyhedral Gravity Dependencies
#    Description: Dependencies required for the polyhedral gravity modelling
#==============================================================================
POLYHEDRAL_GRAVITY_DEPS = [
    {
        "type": git_repository,
        "name": "com_github_nvidia_thrust",
        "remote": "https://github.com/NVIDIA/thrust.git",
        "commit": "3cd56842c94de4926157f6ccdfbbf03ef7e5d5dc",
        "init_submodules": True,
        "build_file": "@//odin:external/thrust.BUILD",
    },
    {
        "type": http_archive,
        "name": "com_github_jbeder_yaml_cpp",
        "urls": ["https://github.com/jbeder/yaml-cpp/archive/refs/tags/yaml-cpp-0.7.0.zip"],
        "strip_prefix": "yaml-cpp-yaml-cpp-0.7.0",
        "sha256": "4d5e664a7fb2d7445fc548cc8c0e1aa7b1a496540eb382d137e2cc263e6d3ef5",
    },
    {
        "type": git_repository,
        "name": "com_github_xtensor_stack_xsimd",
        "remote": "https://github.com/xtensor-stack/xsimd.git",
        "commit": "1577b02d549cca52aa5e943c16f2600950480289",
        "build_file": "@//odin:external/xsimd.BUILD",
    },
    {
        "type": git_repository,
        "name": "com_github_libigl_tetgen",
        "remote": "https://github.com/libigl/tetgen.git",
        "commit": "4f3bfba3997f20aa1f96cfaff604313a8c2c85b6",
        "build_file": "@//odin:external/tetgen.BUILD",
    },
    {
        "type": http_archive,
        "name": "com_github_oneapi_onetbb",
        "sha256": "fcebb93cb9f7e882f62cd351b1c093dbefdcae04b616227dc716b0a5efa9e8ab",
        "strip_prefix": "oneTBB-2021.9.0",
        "build_file": "@//odin:external/tbb.BUILD",
        "urls": ["https://github.com/oneapi-src/oneTBB/archive/v2021.9.0.zip"],
    },
    {
        "type": git_repository,
        "name": "com_github_esa_polyhedral_gravity_model",
        "remote": "https://github.com/geoffreygarrett/polyhedral-gravity-model.git",
        "commit": "2b218d5ed2943a86f5adf0fa9c1399a4a0358594",
        "build_file": "@//odin:external/polyhedral_gravity_model.BUILD",
    },
]

def polyhedral_gravity_dependencies():
    for dep in POLYHEDRAL_GRAVITY_DEPS:
        maybe(
            dep["type"],
            name = dep["name"],
            remote = dep.get("remote"),
            commit = dep.get("commit"),
            init_submodules = dep.get("init_submodules"),
            build_file = dep.get("build_file"),
            urls = dep.get("urls"),
            strip_prefix = dep.get("strip_prefix"),
            sha256 = dep.get("sha256"),
        )

#==============================================================================
#                               Odin Dependencies
#    Description: General dependencies required for the Odin project
#==============================================================================
def odin_dependencies(rules_foreign_cc = True, spdlog = True, fmtlib = True, thrust = True):
    maybe(
        http_archive,
        name = "rules_pkg",
        sha256 = "360c23a88ceaf7f051abc99e2e6048cf7fe5d9af792690576554a88b2013612d",
        strip_prefix = "rules_pkg-0.9.1",
        urls = ["https://github.com/bazelbuild/rules_pkg/archive/0.9.1.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_gabime_spdlog",
        urls = ["https://github.com/gabime/spdlog/archive/v1.11.0.tar.gz"],
        build_file = "@//odin:external/spdlog.BUILD",
        sha256 = "ca5cae8d6cac15dae0ec63b21d6ad3530070650f68076f3a4a862ca293a858bb",
        strip_prefix = "spdlog-1.11.0",
    )

    maybe(
        http_archive,
        name = "rules_cuda",
        sha256 = "dc1f4f704ca56e3d5edd973f98a45f0487d0f28c689d0a57ba236112148b1833",
        strip_prefix = "rules_cuda-v0.1.2",
        urls = ["https://github.com/bazel-contrib/rules_cuda/releases/download/v0.1.2/rules_cuda-v0.1.2.tar.gz"],
    )

    polyhedral_gravity_dependencies()

    maybe(
        git_repository,
        name = "com_github_nvidia_thrust",
        remote = "https://github.com/NVIDIA/thrust.git",
        commit = "3cd56842c94de4926157f6ccdfbbf03ef7e5d5dc",
        init_submodules = True,  # Make sure to pull in any submodules
        build_file = "//:external/thrust.BUILD",
    )

    maybe(
        http_archive,
        name = "com_github_fmtlib_fmt",
        urls = ["https://github.com/fmtlib/fmt/archive/refs/tags/9.1.0.tar.gz"],
        build_file = "@//odin:external/fmt.BUILD",
        sha256 = "5dea48d1fcddc3ec571ce2058e13910a0d4a6bab4cc09a809d8b1dd1c88ae6f2",
        strip_prefix = "fmt-9.1.0",
    )

    maybe(
        git_repository,
        name = "rules_foreign_cc",
        remote = "https://github.com/bazelbuild/rules_foreign_cc",
        commit = "816905a078773405803e86635def78b61d2f782d",
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
        build_file = "@//odin:external/eigen.BUILD",
        sha256 = "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72",
        strip_prefix = "eigen-3.4.0",
        urls = ["https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_autodiff_autodiff",
        build_file = "@//odin:external/bazel_autodiff.BUILD",
        sha256 = "ce87b642272a84f3dcce5463d0f97ef30b15e5a5645898061c06c16a8d61bf22",
        strip_prefix = "autodiff-1.0.3",
        urls = ["https://github.com/autodiff/autodiff/archive/v1.0.3.zip"],
    )

    maybe(
        http_archive,
        name = "com_github_uscilab_cereal",
        build_file = "@//odin:external/cereal.BUILD",
        sha256 = "e72c3fa8fe3d531247773e346e6824a4744cc6472a25cf9b30599cd52146e2ae",
        strip_prefix = "cereal-1.3.2",
        urls = ["https://github.com/USCiLab/cereal/archive/v1.3.2.zip"],
    )

    http_archive(
        name = "build_bazel_apple_support",
        sha256 = "45d6bbad5316c9c300878bf7fffc4ffde13d620484c9184708c917e20b8b63ff",
        url = "https://github.com/bazelbuild/apple_support/releases/download/1.8.1/apple_support.1.8.1.tar.gz",
    )

    maybe(
        git_repository,
        name = "rules_metal",
        remote = "https://github.com/nchavez324/rules_metal",
        commit = "96594931692cef6d842a7a2dbc9249882d838bf3",
    )

    maybe(
        http_archive,
        name = "developer_apple_com_metal_cpp",
        build_file = "@//odin:external/metal_cpp.BUILD",
        urls = ["https://developer.apple.com/metal/cpp/files/metal-cpp_macOS14_iOS17-beta.zip"],
        strip_prefix = "metal-cpp",
        sha256 = "2009a339ecbd56b36601435fe08c415749f8ad09145755472bb637b319003367",
    )

    #    http_archive(
    #        name = "metal_cpp_macOS13_iOS16",
    #        urls = ["https://developer.apple.com/metal/cpp/files/metal-cpp_macOS13_iOS16.zip"],
    #        build_file = "@//:BUILD.metal_cpp_macOS13_iOS16",
    #    )

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

    native.new_local_repository(
        name = "local_brew_gnu_gsl",
        build_file_content = """
cc_library(
    name = "gsl",
    hdrs = glob(["include/**/*.h"]),
    includes = ["include"],
    visibility = ["//visibility:public"],
)
""",
        path = "/opt/homebrew/opt/gsl/",
    )

    native.new_local_repository(
        name = "local_apt_gnu_gsl",
        build_file_content = """
cc_library(
    name = "gsl",
    hdrs = glob(["include/**/*.h"]),
    includes = ["include"],
    visibility = ["//visibility:public"],
)
""",
        path = "/usr/include/gsl/",
    )

    maybe(
        http_archive,
        name = "org_gnu_gsl",
        build_file = "@//odin:external/gsl.BUILD",
        sha256 = "dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b",
        strip_prefix = "gsl-2.7.1",
        urls = ["https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_oneapi_onetbb",
        sha256 = "fcebb93cb9f7e882f62cd351b1c093dbefdcae04b616227dc716b0a5efa9e8ab",
        strip_prefix = "oneTBB-2021.9.0",
        build_file = "@//odin:external/tbb.BUILD",
        urls = ["https://github.com/oneapi-src/oneTBB/archive/v2021.9.0.zip"],
    )

    maybe(
        http_archive,
        name = "com_github_bluebrain_highfive",
        build_file = "@//odin:external/highfive.BUILD",
        sha256 = "029fb6b13fe6ef8098e879e95a55f501f34760d2f9aaa2bdfd513bba1058276b",
        strip_prefix = "HighFive-2.7.1",
        urls = ["https://github.com/BlueBrain/HighFive/archive/v2.7.1.zip"],
    )

    maybe(
        http_archive,
        name = "net_zlib",
        build_file = "@//odin:external/zlib.BUILD",
        sha256 = "b3a24de97a8fdbc835b9833169501030b8977031bcb54b3b3ac13740f846ab30",
        strip_prefix = "zlib-1.2.13",
        urls = ["https://zlib.net/zlib-1.2.13.tar.gz"],
    )

    maybe(
        http_archive,
        name = "com_github_hdfgroup_hdf5",
        build_file = "@//odin:external/hdf5.BUILD",
        sha256 = "5fbebeb7162d336ed06f119ba32f336f2a67c0a492c7802e48f16d4cfb8b0fbf",
        strip_prefix = "hdf5-hdf5-1_14_1-2",
        urls = ["https://github.com/HDFGroup/hdf5/archive/hdf5-1_14_1-2.zip"],
    )
