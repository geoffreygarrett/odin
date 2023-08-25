load("@bazel_tools//tools/build_defs/repo:utils.bzl", "maybe")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def _gsl_repository_impl(repository_ctx):
    gsl_path = repository_ctx.attr.bazel_gsl_path

    if gsl_path:
        native.new_local_repository(
            name = repository_ctx.name,
            build_file_content = """
cc_library(
    name = "gsl",
    hdrs = glob(["include/**/*.h"]),
    includes = ["include"],
    visibility = ["//visibility:public"],
)
""",
            path = gsl_path,
        )
    else:
        maybe(
            http_archive,
            name = repository_ctx.name,
            build_file = "@//odin:external/gsl.BUILD",
            sha256 = "dcb0fbd43048832b757ff9942691a8dd70026d5da0ff85601e52687f6deeb34b",
            strip_prefix = "gsl-2.7.1",
            urls = ["https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz"],
        )

gsl_repository = repository_rule(
    implementation = _gsl_repository_impl,
    attrs = {
        "bazel_gsl_path": attr.string(),
    },
    local = True,
)
