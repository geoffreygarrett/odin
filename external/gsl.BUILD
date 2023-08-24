load("@rules_foreign_cc//foreign_cc:defs.bzl", "configure_make")

#config_setting(
#    name = "linux",
#    values = {"cpu": "k8"},
#    visibility = ["//visibility:public"],
#)
#
#config_setting(
#    name = "macos",
#    values = {"cpu": "darwin"},
#    visibility = ["//visibility:public"],
#)
#
#config_setting(
#    name = "windows",
#    values = {"cpu": "x64_windows"},
#    visibility = ["//visibility:public"],
#)

configure_make(
    name = "gsl",
    lib_source = "//:gsl_sources",
    configure_options = select({
        "@platforms//os:windows": ["--disable-static"],
        "@platforms//os:macos": ["--program-prefix=g", "--disable-static"],
        "//conditions:default": [],
    }),
    #    }),
    copts = select({
        "@platforms//os:windows": ["-DGSL_DLL", "-DWIN32"],
        "//conditions:default": [],
    }),
    args = ["-j$(nproc)"],
    #    linkopts = select({
    #        "@platforms//os:windows": ["-lcblas"],
    #        #        "//conditions:default": ["-lcblas", "-lm"],
    #    }),
    out_include_dir = "include",
    out_lib_dir = "lib",
    out_static_libs = select({
        "@platforms//os:linux": ["libgsl.a", "libgslcblas.a"],
        "@platforms//os:macos": ["libgsl.a", "libgslcblas.a"],
        "@platforms//os:windows": ["gsl.lib", "gslcblas.lib"],
        "//conditions:default": [],
    }),
    out_shared_libs = select({
        "@platforms//os:linux": ["libgsl.so", "libgslcblas.so"],
        "@platforms//os:macos": ["libgsl.dylib", "libgslcblas.dylib"],
        "@platforms//os:windows": ["gsl.dll", "gslcblas.dll"],
        "//conditions:default": [],
    }),
    visibility = ["//visibility:public"],
    env = select({
        "@platforms//os:macos": {
            "AR": "$$(command -v glibtool)",  # brew install libtool
        },
        "//conditions:default": {},
    }),
)

filegroup(
    name = "gsl_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)

#filegroup(
#    name = "gsl_brew_headers",
#    srcs = glob(["/opt/homebrew/include/gsl/**"]),
#)
#
#cc_import(
#    name = "gsl",
#    shared_library = ":gsl",
#    static_library = ":gsl",
#    include =
#    visibility = ["//visibility:public"],
#)
