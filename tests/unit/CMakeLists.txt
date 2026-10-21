# Unit tests: problems with known answers, run by ctest (label "unit").
#
#   ctest --test-dir build -L unit --output-on-failure

# rh_unit_test(NAME [GEOMETRY DIR TOPOLOGY T] SOURCES...): a test program
# linked against librh and the global variables of the programs
# (globals.c), with the geometry.h of DIR (default rhf1d, ONE_D_PLANE).
function(rh_unit_test name)
  cmake_parse_arguments(T "" "GEOMETRY;TOPOLOGY" "" ${ARGN})
  if(NOT T_GEOMETRY)
    set(T_GEOMETRY rhf1d)
    set(T_TOPOLOGY ONE_D_PLANE)
  endif()
  add_executable(${name} ${T_UNPARSED_ARGUMENTS} globals.c)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR}/${T_GEOMETRY}
                             ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_definitions(${name} PRIVATE RH_TOPOLOGY=${T_TOPOLOGY})
  target_link_libraries(${name} PRIVATE rh)
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS unit)
endfunction()

rh_unit_test(test_voigt test_voigt.c)
rh_unit_test(test_linalg test_linalg.c)
rh_unit_test(test_coupling test_coupling.c
  ${PROJECT_SOURCE_DIR}/rhf1d/project.c)
rh_unit_test(test_interp test_interp.c)
rh_unit_test(test_quadrature test_quadrature.c)
rh_unit_test(test_functions test_functions.c)
rh_unit_test(test_formal_1d test_formal_1d.c
  ${PROJECT_SOURCE_DIR}/rhf1d/piecewise_1D.c
  ${PROJECT_SOURCE_DIR}/rhf1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/rhf1d/feautrier.c
  ${PROJECT_SOURCE_DIR}/rhf1d/project.c)
rh_unit_test(test_formal_2d GEOMETRY rhsc2d TOPOLOGY TWO_D_PLANE
  test_formal_2d.c
  ${PROJECT_SOURCE_DIR}/rhsc2d/piecewise_2D.c
  ${PROJECT_SOURCE_DIR}/rhsc2d/bezier_2D.c
  ${PROJECT_SOURCE_DIR}/rhsc2d/geometry.c
  ${PROJECT_SOURCE_DIR}/rhsc2d/solvelong.c
  ${PROJECT_SOURCE_DIR}/rhsc2d/piecestokes_2D.c
  ${PROJECT_SOURCE_DIR}/rhsc2d/project.c)
rh_unit_test(test_formal_3d GEOMETRY rhsc3d TOPOLOGY THREE_D_PLANE
  test_formal_3d.c
  ${PROJECT_SOURCE_DIR}/rhsc3d/piecewise_3D.c
  ${PROJECT_SOURCE_DIR}/rhsc3d/bezier_3D.c
  ${PROJECT_SOURCE_DIR}/rhsc3d/geometry.c
  ${PROJECT_SOURCE_DIR}/rhsc3d/solvelong.c
  ${PROJECT_SOURCE_DIR}/rhsc3d/piecestokes_3D.c
  ${PROJECT_SOURCE_DIR}/rhsc3d/project.c)
rh_unit_test(test_formal_sphere GEOMETRY rhsphere TOPOLOGY SPHERICAL_SYMMETRIC
  test_formal_sphere.c
  ${PROJECT_SOURCE_DIR}/rhsphere/feautrier.c
  ${PROJECT_SOURCE_DIR}/rhsphere/geometry.c
  ${PROJECT_SOURCE_DIR}/rhsphere/project.c)
rh_unit_test(test_formal_stokes test_formal_stokes.c
  ${PROJECT_SOURCE_DIR}/rhf1d/piecewise_1D.c
  ${PROJECT_SOURCE_DIR}/rhf1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/rhf1d/project.c)
rh_unit_test(test_zeeman test_zeeman.c
  ${PROJECT_SOURCE_DIR}/rhf1d/project.c
  ${PROJECT_SOURCE_DIR}/rhf1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/rhf1d/piecewise_1D.c)
