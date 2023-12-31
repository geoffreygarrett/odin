load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

filegroup(
    name = "eigen_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:private"],
    # Private is default, but just to be explicit. Sources would only be
    # part of the cmake target, but not part of the cc_library target.
)

filegroup(
    name = "eigen_headers",
    srcs = glob(["Eigen/**", "unsupported/Eigen/**"]),
    visibility = ["//visibility:public"],
    # We make this public, so we have the option of packaging it
    # with our output.
)

cc_library(
    name = "header_lib",
    hdrs = [":eigen_headers"],
    # We don't need to include the sources here, because this library
    # is a header-only library. CMake is just used for testing etc.
    visibility = ["//visibility:public"],
    includes = ["."],
)

# Even though eigen is header only, libraries like pagmo are easier
# to build with a cmake target. So we create a cmake target that
# performs the installation of the headers and cmake configs.
cmake(
    name = "eigen_cmake",
    lib_source = ":eigen_sources",
    includes = [".", "unsupported"],
    out_include_dir = "",
    out_headers_only = True,  # Flag variable to indicate that the library produces only headers
    install = True,
    linkopts = select({
        "//conditions:default": [""],
        "@platforms//os:macos": [
            "-lomp",
            #                                 "-fopenmp"
        ],
    }),
    copts = select({
        "//conditions:default": [""],
        #        "@platforms//os:macos": ["-fopenmp"],
    }),
    deps = [
        ":header_lib",
    ] + select({
        "@platforms//os:osx": ["@libomp//:omp"],
        "//conditions:default": [],
    }),
    cache_entries = {
        #        "CMAKE_INSTALL_DATADIR": "lib","$$EXT_BUILD_DEPS$$/eigen_cmake/share/eigen3/cmake"
        #        "CMAKE_INSTALL_PREFIX": "$$INSTALLDIR$$",
        "INCLUDE_INSTALL_DIR": "$$INSTALLDIR$$/include",
    },
    visibility = ["//visibility:public"],
)

#filegroup(
#     name = "eigen_sources",
#     srcs = glob(["**"]),
#     visibility = ["//visibility:public"],
# )

cc_library(
    name = "eigen",
    #    hdrs = [":eigen_cmake"],
    hdrs = [":header_lib"],
    # We don't need to include the sources here, because this library
    # is a header-only library. CMake is just used for testing etc.
    visibility = ["//visibility:public"],
    linkopts = select({
        "//conditions:default": [""],
        #        "@platforms//os:macos": ["-lomp"],
    }),
    copts = select({
        "//conditions:default": [""],
        #        "@platforms//os:macos": ["-fopenmp"],
    }),
    includes = select({
        "//conditions:default": ["."],
    }),
)
