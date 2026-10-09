cmake_minimum_required(VERSION 3.18)
project(WaveQLab3D LANGUAGES Fortran)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(MPI REQUIRED COMPONENTS Fortran)
include(CTest)
include(CheckFortranCompilerFlag)
include(CheckFortranSourceCompiles)

# Legacy source relies on eight-byte default REAL. Keep this consistent across
# configurations without replacing user-supplied optimization/debug flags.
set(WQL3D_REAL8_FLAG "" CACHE STRING "Override compiler option(s) for eight-byte default REAL")
if(NOT WQL3D_REAL8_FLAG)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(real8_flags -fdefault-real-8 -fdefault-double-8)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    if(WIN32)
      set(real8_flags /real-size:64)
    else()
      set(real8_flags -r8)
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
    set(real8_flags -r8)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Cray")
    set(real8_flags -s real64)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(LLVMFlang|Flang)$")
    set(real8_flags -fdefault-real-8 -fdefault-double-8)
  else()
    message(FATAL_ERROR "Unsupported Fortran compiler ${CMAKE_Fortran_COMPILER_ID}. Set WQL3D_REAL8_FLAG to its eight-byte REAL option(s).")
  endif()
else()
  separate_arguments(real8_flags NATIVE_COMMAND "${WQL3D_REAL8_FLAG}")
endif()
# Key cached checks by the selected flags so changing an override rechecks it.
string(SHA256 flag_key "${real8_flags}")
check_fortran_compiler_flag("${real8_flags}" "WQL3D_HAS_REAL8_${flag_key}")
if(NOT WQL3D_HAS_REAL8_${flag_key})
  message(FATAL_ERROR "Compiler does not support ${real8_flags}; set WQL3D_REAL8_FLAG for this toolchain.")
endif()

set(saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
string(JOIN " " precision_flags ${real8_flags})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${precision_flags}")
check_fortran_source_compiles("
program precision_check
  integer, parameter :: assert_real64 = 1 / merge(1, 0, storage_size(1.0) == 64)
  integer, parameter :: assert_double64 = 1 / merge(1, 0, storage_size(1.0d0) == 64)
  print *, assert_real64, assert_double64
end program
" "WQL3D_PRECISION_${flag_key}" SRC_EXT F90)
set(CMAKE_REQUIRED_FLAGS "${saved_required_flags}")
if(NOT WQL3D_PRECISION_${flag_key})
  message(FATAL_ERROR "The selected flags must make default REAL and DOUBLE PRECISION eight bytes. Set WQL3D_REAL8_FLAG accordingly.")
endif()

function(wql3d_setup_target target)
  set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/${target}_modules")
  file(MAKE_DIRECTORY "${module_dir}")
  set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
  target_include_directories(${target} PRIVATE "${module_dir}")
  # Intel's Visual Studio integration adds a configuration subdirectory.
  if(CMAKE_GENERATOR MATCHES "Visual Studio" AND CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    target_include_directories(${target} PRIVATE "${module_dir}/$<CONFIG>")
  endif()
  target_compile_options(${target} PRIVATE ${real8_flags})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${target} PRIVATE -Wuninitialized
      "$<$<CONFIG:Debug>:-Wall;-Wextra;-fcheck=bounds;-fbacktrace;-finit-real=nan;-ffpe-trap=invalid,zero,overflow>")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
    if(WIN32)
      target_compile_options(${target} PRIVATE "$<$<CONFIG:Debug>:/check:bounds;/traceback>")
    else()
      target_compile_options(${target} PRIVATE "$<$<CONFIG:Debug>:-check;bounds;-traceback>")
    endif()
  endif()
endfunction()

# Setup the source files and the executable name for the project.
set(WQL3D_SRC
  fault_output.f90
  fault_receiver_output.f90
  iface.f90
  block.f90
  boundary.f90
  common.f90
  diagnostics.f90
  decomposition_safety.f90
  datatypes.f90
  domain.f90
  metrics.f90
  elastic.f90
  plastic.f90
  fields.f90
  grid.f90
  main.f90
  withers_tables.f90
  anelastic_q4_model.f90
  anelastic_q8_model.f90
  anelastic_cq8_b2_model.f90
  anelastic_cq_model.f90
  anelastic_fq_model.f90
  anelastic_fq8_model.f90
  simulation_config.f90
  input_preflight.f90
  material.f90
  anelastic_cq_material.f90
  anelastic_fq_material.f90
  plastic_material.f90
  pml.f90
  mpi3dbasic.f90
  mpi3dcomm.f90
  mpi3dio.f90
  mpi3d_interface.f90
  time_step.f90
  JU_xJU_yJU_z6.f90
  RHS_Interior.f90
  mms.f90
  BoundaryConditions.f90
  CouplingForcing.f90
  Interface_Condition.f90
  initial_stress_condition.f90
  inter_material.f90
  unit_normals.f90
  seismogram.f90
  slice_output.f90
  plane_output.f90
  #hdf5_output.f90
  parallel_write.f90
  moment_tensor.f90
  )

set (PRE_WQL3D_SRC
    common.f90
    diagnostics.f90
    datatypes.f90
    grid.f90
    preprocessor.f90
    metrics.f90
    mpi3dbasic.f90
    mpi3dcomm.f90
    mpi3dio.f90
    mpi3d_interface.f90
  )

add_executable (waveqlab3d ${WQL3D_SRC})
add_executable (pre_wql3d ${PRE_WQL3D_SRC})
foreach(target waveqlab3d pre_wql3d)
  wql3d_setup_target(${target})
  target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
endforeach()
install(TARGETS waveqlab3d pre_wql3d RUNTIME DESTINATION bin)

################################################
# Enable testing using example files
###############################################

if(BUILD_TESTING)

add_executable(fq_coefficients_test
  common.f90 withers_tables.f90 anelastic_fq_model.f90
  anelastic_cq_model.f90 anelastic_fq8_model.f90 ../tests/fq_coefficients_test.f90)
wql3d_setup_target(fq_coefficients_test)
add_test(NAME fq_coefficients_unit COMMAND $<TARGET_FILE:fq_coefficients_test>)

add_executable(fq8_effective_response_test
  common.f90
  withers_tables.f90
  anelastic_fq8_model.f90
  ../tests/fq8_effective_response_test.f90)
wql3d_setup_target(fq8_effective_response_test)
add_test(NAME fq8_effective_response_unit
  COMMAND $<TARGET_FILE:fq8_effective_response_test>)

add_executable(cq_coefficients_test
  common.f90
  anelastic_cq_model.f90
  ../tests/cq_coefficients_test.f90)
wql3d_setup_target(cq_coefficients_test)
add_test(NAME cq_coefficients_unit COMMAND $<TARGET_FILE:cq_coefficients_test>)

# Override for cluster wrappers when the launcher found by FindMPI is unsuitable.
set(WQL3D_MPI_RUNNER "${MPIEXEC_EXECUTABLE}" CACHE FILEPATH "MPI test launcher")
if(WQL3D_MPI_RUNNER AND EXISTS "${WQL3D_MPI_RUNNER}")
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../cmake/mpi_test_config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake" @ONLY)
else()
  message(WARNING "MPI launcher unavailable; only coefficient unit tests will be registered.")
endif()
find_package(Python3 QUIET COMPONENTS Interpreter NumPy)
if(NOT Python3_Interpreter_FOUND OR NOT Python3_NumPy_FOUND)
  message(STATUS "Legacy binary comparison tests require Python 3 and NumPy; skipping them.")
endif()
set(legacy_test_args
  "-DROOT=${CMAKE_CURRENT_SOURCE_DIR}/.."
  "-DEXE=$<TARGET_FILE:waveqlab3d>"
  "-DPRE_EXE=$<TARGET_FILE:pre_wql3d>"
  "-DPYTHON=${Python3_EXECUTABLE}")
if(WQL3D_MPI_RUNNER AND EXISTS "${WQL3D_MPI_RUNNER}")

set (no_mpi_proc 4)

file(MAKE_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/test/1/data
  ${CMAKE_CURRENT_BINARY_DIR}/test/2/data
  ${CMAKE_CURRENT_BINARY_DIR}/test/3/data
  ${CMAKE_CURRENT_BINARY_DIR}/test/4/data)

if(Python3_Interpreter_FOUND AND Python3_NumPy_FOUND)
add_test(NAME mpi_test1
  COMMAND ${CMAKE_COMMAND}
    ${legacy_test_args}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-Din=test_rup_cart_fric.in"
    "-Dprefix=test_rup_cart_fric"
    "-Dt=1"
    "-Dn=${no_mpi_proc}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)
endif()

add_test(NAME fq_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ_dynamic.in"
    "-DSTATE_LABEL=fQ"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME fq_mechanisms_4_through_8
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ_dynamic.in"
    "-DSTATE_LABEL=fQ"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_fq_mechanism_regression.cmake)

add_test(NAME fq_dynamic_pml_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ_dynamic.in"
    "-DSTATE_LABEL=fQ"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_fq_pml_regression.cmake)

add_test(NAME fq_traditional_and_drp_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ_dynamic.in"
    "-DSTATE_LABEL=fQ"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_fq_fd_regression.cmake)

add_test(NAME fq_rejects_invalid_parameters
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ_dynamic.in"
    "-DSTATE_LABEL=fQ"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_fq_invalid_regression.cmake)

add_test(NAME q8_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_dynamic.in"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME cq8_b2_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ8_b2_dynamic.in"
    "-DSTATE_LABEL=cQ8-b2"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME cq8_b2_rejects_one_block
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DNPROCS=1"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ8_b2_invalid_nb1.in"
    "-DEXPECTED_CODE=CFG-CQ8-B2-002"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)

add_test(NAME cq_mechanisms_4_through_8
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_mechanism_regression.cmake)

add_test(NAME cq_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in"
    "-DSTATE_LABEL=cQ"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME cq_rejects_invalid_mechanism_counts
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_invalid_regression.cmake)

add_test(NAME cq_dynamic_pml_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_pml_regression.cmake)

add_test(NAME cq_traditional_and_drp_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_cQ_dynamic.in"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_cq_fd_regression.cmake)

add_test(NAME q4_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_dynamic.in"
    "-DSTATE_LABEL=Q4"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME q4_dynamic_pml_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_dynamic_pml.in"
    "-DSTATE_LABEL=Q4"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME q4_asymmetric_three_ranks
  COMMAND "${WQL3D_MPI_RUNNER}" ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:waveqlab3d> ${MPIEXEC_POSTFLAGS}
    ${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_asymmetric.in)

add_test(NAME q4_deprecated_alias
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_alias.in"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q4_alias_regression.cmake)

foreach(q4_invalid_case legacy partial method nonpositive band)
  add_test(NAME q4_rejects_${q4_invalid_case}
    COMMAND ${CMAKE_COMMAND}
      "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
      "-DNPROCS=2"
      "-DEXE=$<TARGET_FILE:waveqlab3d>"
      "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q4_invalid_${q4_invalid_case}.in"
      "-DEXPECTED_CODE=CFG-Q4-001"
      -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)
endforeach()

add_test(NAME q8_dynamic_pml_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_dynamic_pml.in"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME fq8_coarse_grain_2_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ8_withers.in"
    "-DSTATE_LABEL=fQ8"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME fq8_coarse_grain_0_dynamic_decomposition
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ8_full.in"
    "-DSTATE_LABEL=fQ8"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_q8_dynamic_regression.cmake)

add_test(NAME fq8_rejects_full_withers
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DNPROCS=2"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_fQ8_invalid_full_withers.in"
    "-DEXPECTED_CODE=CFG-FQ8-001"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)

add_test(NAME q8_asymmetric_three_ranks
  COMMAND "${WQL3D_MPI_RUNNER}" ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:waveqlab3d> ${MPIEXEC_POSTFLAGS}
    ${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_asymmetric.in)

add_test(NAME decomposition_rejects_excessive_split
  COMMAND ${CMAKE_COMMAND}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-DNPROCS=4"
    "-DEXE=$<TARGET_FILE:waveqlab3d>"
    "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/../inputfile/test_anelastic_Q8_current.in"
    "-DEXPECTED_CODE=CFG-DECOMP-101"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_expected_preflight_failure.cmake)

if(Python3_Interpreter_FOUND AND Python3_NumPy_FOUND)
add_test(NAME mpi_test2
  COMMAND ${CMAKE_COMMAND}
    ${legacy_test_args}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-Din=test_rup_curv_fric.in"
    "-Dprefix=test_rup_curv_fric"
    "-Dt=2"
    "-Dn=${no_mpi_proc}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)
endif()

if(Python3_Interpreter_FOUND AND Python3_NumPy_FOUND)
add_test(NAME premesh_test1
  COMMAND ${CMAKE_COMMAND}
    ${legacy_test_args}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-Din=test_rup_cart_fric_premesh.in"
    "-Dprefix=test_rup_cart_fric"
    "-Dt=1"
    "-Dn=${no_mpi_proc}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_test.cmake)
endif()

if(Python3_Interpreter_FOUND AND Python3_NumPy_FOUND)
add_test(NAME premesh_test2
  COMMAND ${CMAKE_COMMAND}
    ${legacy_test_args}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-Din=test_rup_curv_fric_premesh.in"
    "-Dprefix=test_rup_curv_fric"
    "-Dt=2"
    "-Dn=${no_mpi_proc}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_test.cmake)
endif()

if(Python3_Interpreter_FOUND AND Python3_NumPy_FOUND)
add_test(NAME serial_test1
  COMMAND ${CMAKE_COMMAND}
    ${legacy_test_args}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-Din=test_rup_cart_fric.in"
    "-Dprefix=test_rup_cart_fric"
    "-Dt=1"
    "-Dn=1"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)
endif()

if(Python3_Interpreter_FOUND AND Python3_NumPy_FOUND)
add_test(NAME serial_test2
  COMMAND ${CMAKE_COMMAND}
    ${legacy_test_args}
    "-DMPI_CONFIG=${CMAKE_CURRENT_BINARY_DIR}/mpi_test_config.cmake"
    "-Din=test_rup_curv_fric.in"
    "-Dprefix=test_rup_curv_fric"
    "-Dt=2"
    "-Dn=1"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/run_mpi_test.cmake)
endif()

endif() # MPI launcher
# Each test gets its own outputs, including when CTest runs in parallel.
get_property(wql3d_tests DIRECTORY PROPERTY TESTS)
foreach(test_name IN LISTS wql3d_tests)
  set(test_dir "${CMAKE_CURRENT_BINARY_DIR}/test_runs/${test_name}")
  file(MAKE_DIRECTORY "${test_dir}/data")
  set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${test_dir}")
endforeach()
endif() # BUILD_TESTING

##############################################################
#
# Add a target to generate API documentation with Doxygen
#
##############################################################

find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../conf/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(docs
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)
