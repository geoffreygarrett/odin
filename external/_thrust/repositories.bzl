load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:utils.bzl", "maybe")

def thrust_dependencies(rules_foreign_cc = True, spdlog = True, fmtlib = True):
    maybe(
        http_archive,
        name = "rules_cuda",
        sha256 = "dc1f4f704ca56e3d5edd973f98a45f0487d0f28c689d0a57ba236112148b1833",
        strip_prefix = "rules_cuda-v0.1.2",
        urls = ["https://github.com/bazel-contrib/rules_cuda/releases/download/v0.1.2/rules_cuda-v0.1.2.tar.gz"],
    )
