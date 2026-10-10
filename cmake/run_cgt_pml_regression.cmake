include("${MPI_CONFIG}")
file(READ "${INPUT}" input_text)
string(REPLACE "btp(1)%pml_lqrs=F,F,F" "btp(1)%pml_lqrs=T,F,F" pml_text "${input_text}")
string(REPLACE "btp(2)%pml_rqrs=F,F,F" "btp(2)%pml_rqrs=T,F,F" pml_text "${pml_text}")
string(REPLACE "%npml=0" "%npml=5" pml_text "${pml_text}")
string(REPLACE "t_final=0.02d0" "t_final=0.1d0" pml_text "${pml_text}")
string(REPLACE "btp(2)%rho_s_p=2.7d0,3.464d0,6d0" "btp(2)%rho_s_p=3d0,3.8d0,6.6d0" pml_text "${pml_text}")
string(REPLACE "btp(1)%lqrs=1,1,1" "btp(1)%lqrs=1,2,1" pml_text "${pml_text}")
set(pml_input "${CMAKE_CURRENT_BINARY_DIR}/test_cgt_dynamic_pml.in")
file(WRITE "${pml_input}" "${pml_text}")
execute_process(COMMAND "${CMAKE_COMMAND}"
  "-DMPI_CONFIG=${MPI_CONFIG}" "-DEXE=${EXE}" "-DINPUT=${pml_input}" "-DSTATE_LABEL=CG-T"
  -P ${CMAKE_CURRENT_LIST_DIR}/run_q8_dynamic_regression.cmake
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "CG-T PML regression failed: ${error}\n${output}")
endif()
