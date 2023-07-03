# @article{Biscani2020,
#   doi = {10.21105/joss.02338},
#   url = {https://doi.org/10.21105/joss.02338},
#   year = {2020},
#   publisher = {The Open Journal},
#   volume = {5},
#   number = {53},
#   pages = {2338},
#   author = {Francesco Biscani and Dario Izzo},
#   title = {A parallel global multiobjective framework for optimization: pagmo},
#   journal = {Journal of Open Source Software}
# }
load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

cmake(
    name = "pagmo_cmake",
    lib_source = "//:pagmo_sources",  # This needs to be updated to include your pagmo sources.
    out_include_dir = "include",
    install = True,
    cache_entries = {
        "CMAKE_BUILD_TYPE": "Release",
        "PAGMO_BUILD_TESTS": "ON",
        "PAGMO_BUILD_TUTORIALS": "ON",
        "PAGMO_WITH_EIGEN3": "ON",
        "PAGMO_WITH_NLOPT": "ON",
        "PAGMO_WITH_IPOPT": "OFF",
        #        "PAGMO_WITH_IPOPT": "ON",
        "PAGMO_ENABLE_IPO": "ON",
        "Boost_NO_BOOST_CMAKE": "ON",
        "Eigen3_DIR": "$$EXT_BUILD_DEPS$$/eigen_cmake/share/eigen3/cmake",
    },
    deps = [
        "@com_github_oneapi_onetbb//:tbb",
        "@com_github_eigen_eigen//:eigen_cmake",
        "@com_github_stevengj_nlopt//:nlopt_cmake",

        #        "@com_github_coin_or_ipopt//:ipopt_configure_make",
        # You will need to add your actual dependencies here.
        #        "@boost//:boost",
        #        "@eigen//:eigen",
        #        "@nlopt//:nlopt",
        #        "@ipopt//:ipopt",
    ],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "pagmo_sources",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)

cc_library(
    name = "pagmo",
    deps = [":pagmo_cmake"],  # since pagmo is heavily cmake-bound :'(
    visibility = ["//visibility:public"],
)
