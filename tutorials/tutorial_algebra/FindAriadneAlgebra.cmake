find_library(ARIADNE_ALGEBRA_LIBRARY NAMES ariadne-algebra)

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(GMP QUIET gmp)
  pkg_check_modules(MPFR QUIET mpfr)
endif()

find_path(GMP_INCLUDE_DIR gmp.h HINTS ${GMP_INCLUDE_DIRS})
find_path(MPFR_INCLUDE_DIR mpfr.h HINTS ${MPFR_INCLUDE_DIRS})
find_library(GMP_LIBRARY NAMES gmp libgmp HINTS ${GMP_LIBRARY_DIRS})
find_library(MPFR_LIBRARY NAMES mpfr libmpfr HINTS ${MPFR_LIBRARY_DIRS})
find_path(ARIADNE_ALGEBRA_INCLUDE_DIR ariadne-algebra.hpp PATH_SUFFIXES ariadne-algebra)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneAlgebra DEFAULT_MSG
  ARIADNE_ALGEBRA_LIBRARY
  ARIADNE_ALGEBRA_INCLUDE_DIR
  GMP_INCLUDE_DIR
  MPFR_INCLUDE_DIR
  GMP_LIBRARY
  MPFR_LIBRARY
)

if(AriadneAlgebra_FOUND)
  get_filename_component(ARIADNE_ALGEBRA_INCLUDE_PARENT_DIR ${ARIADNE_ALGEBRA_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_ALGEBRA_INCLUDE_DIRS
    ${ARIADNE_ALGEBRA_INCLUDE_PARENT_DIR}
    ${ARIADNE_ALGEBRA_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneAlgebra::ariadne-algebra)
    add_library(AriadneAlgebra::ariadne-algebra UNKNOWN IMPORTED)
    set_target_properties(AriadneAlgebra::ariadne-algebra PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_ALGEBRA_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_ALGEBRA_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${MPFR_LIBRARY};${GMP_LIBRARY}"
    )
  endif()
endif()
