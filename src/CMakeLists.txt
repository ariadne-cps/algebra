add_library(ALGEBRA_SRC OBJECT
    vector.cpp
    matrix.cpp
    series.cpp
    polynomial.cpp
    multi_index.cpp
    expansion.cpp
    differential.cpp
    algebra.cpp
)

target_include_directories(ALGEBRA_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/algebra"
)

if(MSVC)
    # Explicit class instantiation asks MSVC about declared-but-undefined optional members.
    target_compile_options(ALGEBRA_SRC PRIVATE /wd4661)
endif()
