load("@rules_pkg//pkg:tar.bzl", "pkg_tar")
load("@rules_pkg//pkg:zip.bzl", "pkg_zip")

cc_binary(
    name = "simulate",
    srcs = ["main.cpp"],
    deps = [
        "//:glog",  #TODO: Refactor logging for the in-house macro
        "//:odin",
    ],
)

#python = use_extension("@rules_python//python:extensions.bzl", "python")
load("@bazel_tools//tools/python:toolchain.bzl", "py_runtime_pair")

#py_runtime(
#    name = "python",
#    interpreter_path = "/home/geoffrey/miniconda3/envs/tudat-bundle/bin/python",
#    python_version = "PY3",
#)
#pip_install(
#    name = "numpy",
#    requirements = ":requirements.txt",
#)

#py_runtime_pair(
#    name = "myruntime_pair",
#    py2_runtime = ":none",
#    py3_runtime = ":myruntime",
#)
#
#py_test(
#    name = "plot",
#    srcs = ["plot.py"],
#    data = [
#        ":python",
#        ":simulate",
#    ],  # a cc_binary which we invoke at run time
#    runtime = ":myruntime_pair",
#)

pkg_zip(
    name = "zip",
    srcs = [
        "README.md",
        ":binaries",
    ],
    strip_prefix = "/",
)
