project( ElmerGUI )
SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake/Modules/")
cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Release )
add_definitions ( -w )

#find_package( VTK REQUIRED )
#find_package ( Qt4 REQUIRED )
#find_package( Qwt REQUIRED )
#find_package(OCE REQUIRED)


include ( ${QT_USE_FILE} )
include_directories (
	${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	${QT_QTCORE_INCLUDE_DIR} ${QT_QTGUI_INCLUDE_DIR} ${QT_QTOPENGL_INCLUDE_DIR} ${QT_QTXML_INCLUDE_DIR} ${QT_QTSCRIPT_INCLUDE_DIR}
	.
	../netgen/libsrc/interface
        ${QWT_INCLUDE_DIR}
        ${VTK_INCLUDE_DIRS}
        ${OCC_INCLUDE_DIR}
        ${PYTHON_INCLUDE_DIRS}
	../PythonQt/src)

set ( EXECUTABLE_OUTPUT_PATH . )
set ( ElmerGUI_HDRS)

FILE(GLOB APP_SRC src/*.cpp)
FILE(GLOB APP_HEADERS src/*.h)

IF(USE_VTK)
  FILE(GLOB VTK_SRC vtkpost/*.cpp)
  FILE(GLOB VTK_HEADERS vtkpost/*.h)
  SET(VTK_LIBS QVTK  vtkCommon  vtkDICOMParser 
              vtkFiltering 
              vtkGenericFiltering 
              vtkGraphics 
              vtkHybrid 
              vtkIO 
              vtkImaging 
              vtkInfovis 
              vtkNetCDF 
              vtkRendering 
              vtkViews 
              vtkVolumeRendering 
              vtkWidgets)
ENDIF()

IF(USE_OpenCASCADE)
  FILE(GLOB OCC_SRC cad/*.cpp)
  FILE(GLOB OCC_HEADERS cad/*.h)
ENDIF()

set (ElmerGUI_SRCS ${APP_SRC} ${VTK_SRC} ${OCC_SRC})

FILE(GLOB UI_FILES forms/*.ui vtkpost/*.ui cad/*.ui)
set(ElmerGUI_UIS ${UI_FILES})

QT4_WRAP_UI(UIS ${ElmerGUI_UIS})

set ( ElmerGUI_RSCS ElmerGUI.qrc ElmerGUI.rc)
QT4_ADD_RESOURCES(RSCS ${ElmerGUI_RSCS})

set ( ElmerGUI_MOCS ${APP_HEADERS} ${VTK_HEADERS} ${OCC_HEADERS})
QT4_WRAP_CPP(MOCS ${ElmerGUI_MOCS})

##
SET(ELMERGUI_LINK_LIBS ${ELMERGUI_LINK_LIBS} ${VTK_LIBS} ${PYTHON_LIBRARIES}
  ${OCC_LIBS_TODO} ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY}
  ${QT_QTOPENGL_LIBRARY} ${QT_QTXML_LIBRARY} ${QT_QTSCRIPT_LIBRARY} ${OPENGL_gl_LIBRARY}
  ${OPENGL_glu_LIBRARY})


add_library ( elmergui SHARED ${ElmerGUI_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
target_link_libraries ( elmergui ${ELMERGUI_LINK_LIBS})
ADD_EXECUTABLE(ElmerGUI ${ElmerGUI_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
target_link_libraries ( ElmerGUI ${ELMERGUI_LINK_LIBS})
