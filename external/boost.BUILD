load("@rules_foreign_cc//tools/build_defs:boost_build.bzl", "boost_build")

boost_build(
    name = "boost_fiber",
    lib_source = "@boost//:all",
    static_libraries = ["libboost_fiber.a"],
    user_options = ["--with-fiber"],
    visibility = ["//visibility:public"],
    deps = [":boost_context"],
)

boost_build(
    name = "boost_context",
    lib_source = "@boost//:all",
    static_libraries = ["libboost_context.a"],
    user_options = ["--with-context"],
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "helloworld",
    srcs = [
        "main.cc",
    ],
    deps = [
        ":boost_fiber",
    ],
)
