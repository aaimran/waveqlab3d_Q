include("${MPI_CONFIG}")
set(dir "${CMAKE_CURRENT_BINARY_DIR}/test/${t}")
file(MAKE_DIRECTORY "${dir}/data")
configure_file("${ROOT}/test_problems/${in}" "${dir}/${in}" COPYONLY)
file(COPY "${ROOT}/test_problems/truth/${t}/" DESTINATION "${dir}/truth")
execute_process(COMMAND "${MPIEXEC}" ${MPIEXEC_NUMPROC_FLAG} ${n} ${MPIEXEC_PREFLAGS}
  "${EXE}" ${MPIEXEC_POSTFLAGS} "${in}"
  WORKING_DIRECTORY "${dir}" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Solver failed: ${result}")
endif()
execute_process(COMMAND "${PYTHON}" "${ROOT}/python/read_binary.py" "${dir}" "${prefix}"
  RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "Output comparison failed for ${in}: ${result}")
endif()
