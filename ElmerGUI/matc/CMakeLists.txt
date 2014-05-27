project( matc )
cmake_minimum_required( VERSION 2.6 )
set ( CMAKE_BUILD_TYPE Release )
add_definitions ( -w )
set ( LIBRARY_OUTPUT_PATH lib )
set ( matc_HDRS
	
	)

      FILE(GLOB matc_src_files src/*.c)
set ( matc_SRCS ${matc_src_files})

add_library ( matc SHARED ${matc_SRCS} ${UIS} ${RSCS} ${TRS} ${MOCS} )
target_link_libraries ( matc  )
