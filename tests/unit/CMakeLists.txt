# Unit tests: problems with known answers, run by ctest (label "unit").
#
#   ctest --test-dir build -L unit --output-on-failure

# rh_unit_test(NAME SOURCES...): a test program linked against librh and
# the global variables of the 1-D programs (globals.c).
function(rh_unit_test name)
  add_executable(${name} ${ARGN} globals.c)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR}/rhf1d
                             ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${name} PRIVATE rh)
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS unit)
endfunction()

rh_unit_test(test_voigt test_voigt.c)
rh_unit_test(test_linalg test_linalg.c)
rh_unit_test(test_interp test_interp.c)
rh_unit_test(test_quadrature test_quadrature.c)
rh_unit_test(test_functions test_functions.c)
rh_unit_test(test_formal_1d test_formal_1d.c
  ${PROJECT_SOURCE_DIR}/rhf1d/piecewise_1D.c
  ${PROJECT_SOURCE_DIR}/rhf1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/rhf1d/feautrier.c
  ${PROJECT_SOURCE_DIR}/rhf1d/project.c)
