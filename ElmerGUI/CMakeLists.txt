cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Release )
add_definitions ( -Wall )
SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake/Modules")

find_package ( Qt4 REQUIRED )

OPTION(USE_VTK "Use VTK for postprocessing" OFF)
OPTION(USE_Qwt "Use Qwt for convergence monitor" OFF)
OPTION(USE_OpenCASCADE "Use OpenCASCADE to import cad files" OFF)
OPTION(USE_PythonQT "Use PythonQT for scripting in post processor" OFF)
OPTION(USE_MATC "Use matc for internal operations in postprocessing" OFF)

FIND_PACKAGE(OpenGL REQUIRED)

IF(USE_VTK)
  FIND_PACKAGE(VTK REQUIRED)
  ADD_DEFINITIONS(-DEG_VTKP)
ENDIF()

IF(USE_Qwt)
  FIND_PACKAGE(Qwt REQUIRED)
  ADD_DEFINITIONS(-DEG_QWT)
ENDIF()

IF(USE_OpenCASCADE)
  FIND_PACKAGE(OpenCASCADE REQUIRED)
  ADD_DEFINITIONS(-DEG_OCC)
ENDIF()

IF(USE_PythonQT)
  find_package(PythonLibs REQUIRED)
  ADD_DEFINITIONS(-DEG_PYTHONQT)
ENDIF()

include ( ${QT_USE_FILE} )
include_directories (
  ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  ${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR}
  )

if(USE_MATC)
  add_subdirectory( matc )
endif()

add_subdirectory( netgen )
add_subdirectory( Application )
