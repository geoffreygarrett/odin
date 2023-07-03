## extensions.bzl
load("//:deps.bzl", "odin_dependencies")

def _odin_deps_impl(ctx):
    odin_dependencies()

#
#    for mod in ctx.modules:
#        for dependencies in mod.tags.dependencies:
#            odin_dependencies()

def _dependencies():
    return {}

odin_deps = module_extension(
    #    tag_classes = {
    #        "deps": tag_class(attrs = _dependencies()),
    #    },
    implementation = _odin_deps_impl,
)
