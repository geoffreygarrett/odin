load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

def _impl(ctx):
    print("Starting Miniconda environment setup...")

    # Download and verify Miniconda installer
    print("Downloading and verifying Miniconda installer...")
    result = ctx.execute(["bash", "-c", """
        url="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        hash="aef279d6baea7f67940f16aad17ebe5f6aac97487c7c03466ff01f4819e5a651"
        curl -o miniconda.sh -L $url
        echo "$hash  miniconda.sh" | sha256sum --check
    """])
    print("Miniconda download stdout: ", result.stdout)
    print("Miniconda download stderr: ", result.stderr)
    print("Miniconda download return_code: ", result.return_code)

    # Run Miniconda installer
    print("Running Miniconda installer...")
    result = ctx.execute(["bash", "miniconda.sh", "-b", "-p", "miniconda"])
    print("Miniconda installation stdout: ", result.stdout)
    print("Miniconda installation stderr: ", result.stderr)
    print("Miniconda installation return_code: ", result.return_code)

    # Check if environment file is provided
    if ctx.attr.environment_file:
        print("Installing packages from environment file to base environment...")
        result = ctx.execute(["bash", "-c", """
            export PATH=$PWD/miniconda/bin:$PATH
            ./miniconda/bin/conda env update -n base --file {}
        """.format(ctx.path(ctx.attr.environment_file))])
        print("Environment creation from file stdout: ", result.stdout)
        print("Environment creation from file stderr: ", result.stderr)
        print("Environment creation from file return_code: ", result.return_code)
    else:
        # Or with explicit packages
        print("Installing explicit packages to base environment...")
        package_command = " ".join(ctx.attr.packages)
        result = ctx.execute(["bash", "-c", """
            export PATH=$PWD/miniconda/bin:$PATH
            ./miniconda/bin/conda install -y {}
        """.format(package_command)])
        print("Environment creation with packages stdout: ", result.stdout)
        print("Environment creation with packages stderr: ", result.stderr)
        print("Environment creation with packages return_code: ", result.return_code)

    # Override the dummy BUILD file content with the provided BUILD file
    print("Overriding the dummy BUILD file content with the provided BUILD file...")
    ctx.file("BUILD", ctx.read(ctx.attr.build_file))

    print("Miniconda environment setup completed.")

conda_env = repository_rule(
    implementation = _impl,
    attrs = {
        "env_name": attr.string(),
        "packages": attr.string_list(),
        "environment_file": attr.label(allow_single_file = True),
        "build_file": attr.label(allow_single_file = [".BUILD"]),  # Updated from a list of files to a single file
    },
)

def conda_environment(name, env_name = "", packages = [], environment_file = None, build_file = None):
    conda_env(
        name = name,
        env_name = env_name,
        packages = packages,
        environment_file = environment_file,
        build_file = build_file,  # modified attribute to specify a single BUILD file
    )
