include(${TEST_SOURCE}/../test_macros.cmake)
execute_process(COMMAND ${ELMERGRID_BIN} 14 2 BeamFlow.msh -autoclean )
RUN_ELMER_TEST()
