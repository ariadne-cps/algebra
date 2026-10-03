cmake_minimum_required(VERSION 3.15)

project(Algebra VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/submodules/interval/submodules/numeric/submodules/configuration/cmake/ProjectOptions.cmake")
    message(FATAL_ERROR
        "Git submodules are not initialized. "
        "Run 'git submodule update --init --recursive' before configuring the project.")
endif()

include("${PROJECT_SOURCE_DIR}/submodules/interval/submodules/numeric/submodules/configuration/cmake/ProjectOptions.cmake")

require_initialized_submodules()
setup_project()

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/include)

if(NOT TARGET interval)
    add_subdirectory(submodules)
endif()

if(MSVC AND TARGET pyariadne-numeric-submodule-obj)
    # The standalone pyariadne module embeds NUMERIC_SRC directly, so Numeric
    # Python bindings must reference Numeric symbols as defined by this DLL,
    # not as dllimport symbols from a separate ariadne-numeric DLL.
    target_compile_definitions(pyariadne-numeric-submodule-obj PRIVATE ARIADNE_NUMERIC_BUILD)
endif()

if(MSVC)
    # Numeric headers instantiate conversion paths that MSVC reports as C4244.
    # Keep /WX enabled for Algebra while avoiding warnings owned by the lower layer.
    add_compile_options(/wd4244)
endif()

add_library(algebra INTERFACE)
target_include_directories(algebra INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-algebra>
)
target_link_libraries(algebra INTERFACE interval)

add_subdirectory(src)
target_link_libraries(ALGEBRA_SRC PRIVATE algebra)

ariadne_register_public_headers(
    TARGET algebra
    DIRECTORIES "${PROJECT_SOURCE_DIR}/include/algebra"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-algebra.hpp"
)

setup_project_library(
    TARGET ariadne-algebra
    OBJECTS
        ALGEBRA_SRC
        INTERVAL_SRC
        NUMERIC_SRC
        FOUNDATION_SRC
        UTILITY_SRC
    LINK_LIBRARIES algebra
)

setup_standalone_project_tests(
    TARGET ariadne-algebra
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET algebra
        DESTINATION include/ariadne-algebra
    )
endif()
