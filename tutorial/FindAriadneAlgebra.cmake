find_library(ARIADNE_ALGEBRA_LIBRARY NAMES ariadne-algebra)
find_path(ARIADNE_ALGEBRA_INCLUDE_DIR ariadne-algebra.hpp PATH_SUFFIXES ariadne-algebra)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneAlgebra DEFAULT_MSG
  ARIADNE_ALGEBRA_LIBRARY
  ARIADNE_ALGEBRA_INCLUDE_DIR
)

if(AriadneAlgebra_FOUND)
  get_filename_component(ARIADNE_ALGEBRA_INCLUDE_PARENT_DIR ${ARIADNE_ALGEBRA_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_ALGEBRA_INCLUDE_DIRS
    ${ARIADNE_ALGEBRA_INCLUDE_PARENT_DIR}
    ${ARIADNE_ALGEBRA_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneAlgebra::ariadne-algebra)
    add_library(AriadneAlgebra::ariadne-algebra UNKNOWN IMPORTED)
    set_target_properties(AriadneAlgebra::ariadne-algebra PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_ALGEBRA_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_ALGEBRA_INCLUDE_DIRS}"
    )
  endif()
endif()
