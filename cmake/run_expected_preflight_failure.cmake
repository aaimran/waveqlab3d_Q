include("${MPI_CONFIG}")
execute_process(
  COMMAND "${MPIEXEC}" ${MPIEXEC_NUMPROC_FLAG} "${NPROCS}" ${MPIEXEC_PREFLAGS} "${EXE}" ${MPIEXEC_POSTFLAGS} "${INPUT}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error)
if(result STREQUAL "0")
  message(FATAL_ERROR "invalid input unexpectedly succeeded")
endif()
string(CONCAT combined "${output}" "${error}")
if(NOT combined MATCHES "${EXPECTED_CODE}")
  message(FATAL_ERROR "expected diagnostic ${EXPECTED_CODE} was not emitted: ${combined}")
endif()
