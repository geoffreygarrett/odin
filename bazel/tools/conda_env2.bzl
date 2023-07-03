load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

_OS_INFO = {
    "@bazel_tools//tools/cpp:linux": {
        "url": "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
        "hash": "aef279d6baea7f67940f16aad17ebe5f6aac97487c7c03466ff01f4819e5a651",
    },
    "@bazel_tools//tools/cpp:windows": {
        "url": "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe",
        "hash": "307194e1f12bbeb52b083634e89cc67db4f7980bd542254b43d3309eaf7cb358",
    },
    "@bazel_tools//tools/cpp:darwin": {
        "url": "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh",
        "hash": "5abc78b664b7da9d14ade330534cc98283bb838c6b10ad9cfd8b9cc4153f8104",
    },
    "//conditions:default": {
        "url": "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
        "hash": "aef279d6baea7f67940f16aad17ebe5f6aac97487c7c03466ff01f4819e5a651",
    },
}

def _impl(ctx):
    os_info = select(_OS_INFO)

    # download Miniconda installer
    ctx.download(
        url = os_info["url"],
        sha256 = os_info["hash"],
        output = "miniconda.sh",
    )

    # run Miniconda installer
    ctx.execute(["bash", "miniconda.sh", "-b", "-p", "miniconda"])

    # check if environment file is provided
    if ctx.attr.environment_file:
        ctx.execute(["./miniconda/bin/conda", "env", "create", "-f", ctx.path(ctx.attr.environment_file)])
    else:
        ctx.execute(["./miniconda/bin/conda", "create", "-n", ctx.attr.env_name] + ctx.attr.packages)

conda_env = repository_rule(
    implementation = _impl,
    attrs = {
        "env_name": attr.string(),
        "packages": attr.string_list(),
        "environment_file": attr.label(allow_single_file = True),
    },
)

def conda_environment(name, env_name = "", packages = [], environment_file = None):
    conda_env(
        name = name,
        env_name = env_name,
        packages = packages,
        environment_file = environment_file,
    )
