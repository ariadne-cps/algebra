include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-algebra-submodule-obj OBJECT
    algebra_submodule.cpp
    linear_algebra_submodule.cpp
    polynomial_submodule.cpp
    differentiation_submodule.cpp
)
target_include_directories(pyariadne-algebra-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-interval,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-algebra-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    algebra
)
set_target_properties(pyariadne-algebra-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(pyariadne-algebra INTERFACE)
target_sources(pyariadne-algebra INTERFACE
    $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
)
target_include_directories(pyariadne-algebra INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-interval,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-algebra INTERFACE
    pyariadne-interval
)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    add_library(pyariadne-algebra-module-obj OBJECT
        algebra_module.cpp
    )
    target_include_directories(pyariadne-algebra-module-obj PRIVATE
        "${PYTHON_COMMON_SOURCE_DIR}/include"
    )
    target_link_libraries(pyariadne-algebra-module-obj PRIVATE
        pyariadne-python-common
        pybind11::module
    )
    set_target_properties(pyariadne-algebra-module-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
