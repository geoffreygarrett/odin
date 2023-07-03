#filegroup(
#    name = "include",
#    srcs = glob(["miniconda/envs/odin-env/include/**"]),
#    visibility = ["//visibility:public"],
#)
#
#filegroup(
#    name = "lib",
#    srcs = glob(["miniconda/envs/odin-env/lib/**"]),
#    visibility = ["//visibility:public"],
#)
#
#filegroup(
#    name = "shared",
#    srcs = glob(["miniconda/envs/odin-env/shared/**"]),
#    visibility = ["//visibility:public"],
#)
#filegroup(
#    name = "prefix",
#    srcs = glob(["miniconda/envs/odin-env/**"]),
#    visibility = ["//visibility:public"],
#)

cc_import(
    name = "pagmo",
    hdrs = glob([
        "miniconda/include/pagmo/**/*.hpp",
        "miniconda/include/pagmo/**/*.h",
    ]),
    shared_library = "miniconda/lib/libpagmo.so",
    visibility = ["//visibility:public"],
)

cc_import(
    name = "pagmo",
    hdrs = glob([
        "miniconda/include/pagmo/**/*.hpp",
        "miniconda/include/pagmo/**/*.h",
    ]),
    shared_library = "miniconda/lib/libpagmo.so",
    visibility = ["//visibility:public"],
)

#cc_import(
#    name="prefix",
#
#)

#cc_library(
#    name = "pagmo_cmake",
#    srcs = [],
#    deps = [],
#    visibility = ["//visibility:public"],
#)

#cc_library(
#    name = "boost",
#    hdrs = glob([
#        "miniconda/envs/odin-env/include/boost/**/*.hpp",
#        "miniconda/envs/odin-env/include/boost/**/*.h",
#    ]),
#    srcs = glob([
#        "miniconda/envs/odin-env/**/*/libboost*.so*",
#        "miniconda/envs/odin-env/**/*/libboost*.a",
#    ]),
#    includes = ["miniconda/envs/odin-env/include"],
#    visibility = ["//visibility:public"],
#    linkopts = ["-Wl,-rpath,miniconda/envs/odin-env/lib"],
#)
#
#cc_library(
#    name = "tbb",
#    hdrs = glob([
#        "miniconda/envs/odin-env/include/tbb/**/*.hpp",
#        "miniconda/envs/odin-env/include/tbb/**/*.h",
#    ]),
#    srcs = glob([
#        "miniconda/envs/odin-env/**/*/libtbb*.so*",
#        "miniconda/envs/odin-env/**/*/libtbb*.a",
#    ]),
#    includes = ["miniconda/envs/odin-env/include"],
#    visibility = ["//visibility:public"],
#    linkopts = ["-Wl,-rpath,miniconda/envs/odin-env/lib"],
#)

#def _conda_setup_impl(ctx):
#    pass
#
#conda_setup = rule(
#    implementation = _conda_setup_impl,
#    attrs = {
#        "_conda_env": attr.label(default = Label("@odin_conda_environment//"))
#    },
#)
#
#def conda_environment_setup():
#    conda_setup(name = "conda_setup")

#cc_library(
#    name = "pagmo",
#    hdrs = glob([
#        "miniconda/envs/odin-env/include/pagmo/**/*.hpp",
#        "miniconda/envs/odin-env/include/pagmo/**/*.h",
#        "miniconda/envs/odin-env/include/**/*.h",
#        "miniconda/envs/odin-env/include/**/*.hpp",
#    ]),
#    srcs = glob([
#        "miniconda/envs/odin-env/**/*/libpagmo*.so*",
#        "miniconda/envs/odin-env/**/*.so*",
#        "miniconda/envs/odin-env/**/*.a",
#    ]),
#    #    deps = [":boost", ":tbb"],  # Pagmo depends on Boost
#    includes = ["miniconda/envs/odin-env/include"],
#    visibility = ["//visibility:public"],
#    linkopts = ["-Wl,-rpath,miniconda/envs/odin-env/lib"],
#)
