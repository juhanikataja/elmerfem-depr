INCLUDE(${TEST_SOURCE}/../test_macros.cmake)

EXECUTE_PROCESS(COMMAND ${ELMERGRID_BIN} 1 2 rectangle.grd)

RUN_ELMERICE_TEST()
