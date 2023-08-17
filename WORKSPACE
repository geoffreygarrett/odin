# ---------------------------------------------------------
# Workspace name
# ---------------------------------------------------------
workspace(name = "odin")

#local_repository(
#    name = "pybind11_bazel",
#    path = "external/pybind11_bazel",
#)
#
#local_repository(
#    name = "odin",
#    path = "../odin",
#)

#http_archive(
#    name = "rules_cuda",
#    sha256 = "dc1f4f704ca56e3d5edd973f98a45f0487d0f28c689d0a57ba236112148b1833",
#    strip_prefix = "rules_cuda-v0.1.2",
#    urls = ["https://github.com/bazel-contrib/rules_cuda/releases/download/v0.1.2/rules_cuda-v0.1.2.tar.gz"],
#)
load("@//:repositories.bzl", "odin_dependencies")

odin_dependencies()

load("@rules_cuda//cuda:repositories.bzl", "register_detected_cuda_toolchains", "rules_cuda_dependencies")

rules_cuda_dependencies()

register_detected_cuda_toolchains()

#local_repository(
#    name = "com_github_esa_polyhedral_gravity_model",
#    path = "external/polyhedral-gravity-model",
#)

#alias(
#    name = "com_github_esa_polyhedral_gravity_model",
#    actual = "//:external/polyhedral-gravity-model",
#)

#local_repository(
#    name = "pydin",
#    path = ".",
#)

#load("@pydin//:repositories.bzl", "pydin_dependencies")

#pydin_dependencies()

load("@rules_foreign_cc//foreign_cc:repositories.bzl", "rules_foreign_cc_dependencies")
load("@rules_pkg//pkg:deps.bzl", "rules_pkg_dependencies")

rules_foreign_cc_dependencies()

rules_pkg_dependencies()

#register_execution_platforms(
#    ":x64_windows-clang-cl",
#)
#
#register_toolchains(
#    "@local_config_cc//:cc-toolchain-x64_windows-clang-cl",
#)

#load("@rules_python//python:repositories.bzl", "python_register_toolchains")
#
#python_register_toolchains(
#    name = "python3_10",
#    python_version = "3.10",
#)
#
#load("@rules_python//python:pip.bzl", "pip_install", "pip_parse")
#load("@pybind11_bazel//:python_configure.bzl", "python_configure")
#load("@python3_10//:defs.bzl", "interpreter")
#load("@rules_python//python:pip.bzl", "pip_install", "pip_parse")
#
#python_configure(
#    name = "local_config_python",
#    python_interpreter_target = interpreter,
#)
#
#load("@rules_python//python:pip.bzl", "pip_parse")
#
#pip_parse(
#    name = "pip",
#
#    # Here, we use the interpreter constant that resolves to the host interpreter from the default Python toolchain.
#    python_interpreter_target = interpreter,
#
#    # Uses the default repository name "pip"
#    requirements_lock = "//:requirements_lock.txt",
#    requirements_windows = "//:requirements_windows.txt",
#)
#
#load("@pip//:requirements.bzl", "install_deps")
#
## Initialize repositories for all packages in requirements.txt.
#install_deps()
