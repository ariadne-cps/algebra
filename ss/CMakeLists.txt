if(NOT TARGET interval)
    ariadne_add_submodule(interval)
endif()

# Testing is a test-only dependency, never part of the production library.
ariadne_setup_testing()
