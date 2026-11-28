set(BENCHMARKS benchmark_eigen_decomposition)
if(NOT TARGET benchmarks)
    add_custom_target(benchmarks)
endif()

# Standalone Algebra bundles its functional dependencies, not optional Logging.
# Follow Numeric: inject only Logging objects absent from the Algebra bundle.
get_target_property(_benchmark_logging_objects logging ARIADNE_PROJECT_OBJECT_TARGETS)
get_target_property(_algebra_objects algebra ARIADNE_PROJECT_OBJECT_TARGETS)
if(NOT _benchmark_logging_objects OR NOT _algebra_objects)
    message(FATAL_ERROR "Benchmarking requires registered Algebra/Logging objects")
endif()
foreach(_object IN LISTS _algebra_objects)
    list(REMOVE_ITEM _benchmark_logging_objects "${_object}")
endforeach()
foreach(BENCHMARK IN LISTS BENCHMARKS)
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE
        ariadne-algebra algebra ariadne-benchmarking)
    foreach(_object IN LISTS _benchmark_logging_objects)
        target_sources(${BENCHMARK} PRIVATE $<TARGET_OBJECTS:${_object}>)
    endforeach()
endforeach()
add_dependencies(benchmarks ${BENCHMARKS})
