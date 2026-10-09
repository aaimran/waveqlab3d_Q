include("${MPI_CONFIG}")
file(READ "${INPUT}" input_text)
foreach(fd_type traditional upwind_drp)
  string(REPLACE "fd_type='upwind'" "fd_type='${fd_type}'" case_text "${input_text}")
  set(case_input "${CMAKE_CURRENT_BINARY_DIR}/test_anelastic_fQ_${fd_type}.in")
  file(WRITE "${case_input}" "${case_text}")
  execute_process(COMMAND "${CMAKE_COMMAND}"
    "-DMPI_CONFIG=${MPI_CONFIG}" "-DEXE=${EXE}" "-DINPUT=${case_input}" "-DSTATE_LABEL=fQ"
    -P ${CMAKE_CURRENT_LIST_DIR}/run_q8_dynamic_regression.cmake
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT result STREQUAL "0")
    message(FATAL_ERROR "anelastic-fQ ${fd_type} regression failed: ${error}\n${output}")
  endif()
endforeach()
