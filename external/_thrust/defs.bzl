# Defines
DEFINES = [
    "THRUST_HOST_SYSTEM=CPP",
    "THRUST_DEVICE_SYSTEM=TBB",
]

# Cache Entries
CACHE_ENTRIES = {
    "THRUST_HOST_SYSTEM": "CPP",
    "THRUST_DEVICE_SYSTEM": "TBB",
}

# Exporting variables
def get_cmake_defines():
    return DEFINES

def get_cmake_cache_entries():
    return CACHE_ENTRIES
