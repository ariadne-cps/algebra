add_library(ALGEBRA_SRC OBJECT
    vector.cpp
    matrix.cpp
    series.cpp
    polynomial.cpp
    chebyshev_polynomial.cpp
    multi_index.cpp
    expansion.cpp
    sweeper.cpp
    differential.cpp
    algebra.cpp
)

target_include_directories(ALGEBRA_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/algebra"
)

