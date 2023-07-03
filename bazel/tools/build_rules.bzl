# build_rules.bzl
def create_alias(name, actual, visibility = ["//visibility:public"]):
    native.alias(
        name = name,
        actual = actual,
        visibility = visibility,
    )

def create_aliases(mapping, visibility = ["//visibility:public"]):
    for name, actual in mapping.items():
        create_alias(name, actual, visibility)

def create_config_setting(name, value, visibility = ["//visibility:public"]):
    native.config_setting(
        name = name,
        values = {"define": value},
        visibility = visibility,
    )
