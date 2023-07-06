## extensions.bzl
load("//:repositories.bzl", "odin_dependencies")

def _init(ctx):
    odin_dependencies(rules_foreign_cc = False)

init = module_extension(
    implementation = _init,
)
