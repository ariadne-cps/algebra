if(NOT TARGET interval)
    ensure_submodule_initialized(interval)
    add_subdirectory(interval)
endif()
