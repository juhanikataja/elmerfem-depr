
# CMake script to find pcre library
# Juhani Kataja, CSC - IT Center for Science Ltd.
# 2014/08
#
# Hint variables
# * PCRE_ROOT (env, cmake)
# * PCRE_INCLUDE_DIR (cmake)
# * PCRE_LIBRARY_DIR (cmake)
#
# Create variables PCRE_LIBRARIES and PCRE_INCLUDE_DIR

if(PCRE_LIBRARY)
  if(PCRE_INCLUDE_DIR)
    set(PCRE_FOUND TRUE)
    return()
  endif()
endif()

set(PCRE_FOUND FALSE)

find_path(PCRE_INCLUDE_DIR NAMES pcre.h
  HINTS
  "${PCRE_ROOT}/include"
  "$ENV{PCRE_ROOT}/include"
  ${PCRE_INCLUDE_DIR}
  )

find_library(PCRE_LIBRARY NAMES pcre
  HINTS
  "${PCRE_ROOT}/lib"
  "$ENV{PCRE_ROOT}/lib" 
  "${PCRE_LIBRARY_DIR}")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PCRE DEFAULT_MSG PCRE_LIBRARY PCRE_INCLUDE_DIR)

SET(PCRE_LIBRARIES ${PCRE_LIBRARY})
message(STATUS "PCRE_LIBRARY: ${PCRE_LIBRARY}")
message(STATUS "PCRE_LIBRARIES: ${PCRE_LIBRARIES}")
MARK_AS_ADVANCED(PCRE_INCLUDE_DIR PCRE_LIBRARY PCRE_LIBRARIES)
