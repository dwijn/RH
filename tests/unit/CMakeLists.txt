# Unit tests: problems with known answers, run by ctest (label "unit").
#
#   ctest --test-dir build -L unit --output-on-failure

# rh_unit_test(NAME [GEOMETRY DIR TOPOLOGY T] SOURCES...): a test program
# linked against librh and the global variables of the programs
# (globals.c), with the geometry.h of src/geometry/DIR (default plane1d,
# ONE_D_PLANE).
function(rh_unit_test name)
  cmake_parse_arguments(T "" "GEOMETRY;TOPOLOGY" "" ${ARGN})
  if(NOT T_GEOMETRY)
    set(T_GEOMETRY plane1d)
    set(T_TOPOLOGY ONE_D_PLANE)
  endif()
  add_executable(${name} ${T_UNPARSED_ARGUMENTS} globals.c)
  target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR}/src/geometry/${T_GEOMETRY}
                             ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_definitions(${name} PRIVATE RH_TOPOLOGY=${T_TOPOLOGY})
  target_link_libraries(${name} PRIVATE rh)
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS unit)
endfunction()

rh_unit_test(test_voigt test_voigt.c)
rh_unit_test(test_linalg test_linalg.c)
rh_unit_test(test_gth test_gth.c)
rh_unit_test(test_ng test_ng.c)
rh_unit_test(test_anderson test_anderson.c)
rh_unit_test(test_mpe test_mpe.c)
rh_unit_test(test_negtau test_negtau.c)
rh_unit_test(test_blends test_blends.c ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_mirror test_mirror.c ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_negpop test_negpop.c)
rh_unit_test(test_collision_tables test_collision_tables.c)
rh_unit_test(test_rlklabel test_rlklabel.c ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_prdh test_prdh.c)
rh_unit_test(test_jgas_split test_jgas_split.c)
rh_unit_test(test_newton test_newton.c)
rh_unit_test(test_coupling test_coupling.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_interp test_interp.c)
rh_unit_test(test_quadrature test_quadrature.c)
rh_unit_test(test_stokes_step test_stokes_step.c)
rh_unit_test(test_functions test_functions.c)
rh_unit_test(test_prdfold test_prdfold.c)
rh_unit_test(test_prdrho test_prdrho.c)
rh_unit_test(test_collision_formulas test_collision_formulas.c)
rh_unit_test(test_background test_background.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_hminus_bf test_hminus_bf.c)
rh_unit_test(test_chemequil test_chemequil.c)
target_compile_definitions(test_hminus_bf PRIVATE
  RH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
rh_unit_test(test_hminus_ff test_hminus_ff.c)
rh_unit_test(test_rayleigh_h2 test_rayleigh_h2.c)
rh_unit_test(test_h2ff test_h2ff.c)
target_compile_definitions(test_h2ff PRIVATE
  RH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
rh_unit_test(test_gaunt_ff test_gaunt_ff.c)
target_compile_definitions(test_gaunt_ff PRIVATE
  RH_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
rh_unit_test(test_formal_1d test_formal_1d.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/piecewise_1D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/feautrier.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_bgscatter test_bgscatter.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/piecewise_1D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_formal_2d GEOMETRY cart2d TOPOLOGY TWO_D_PLANE
  test_formal_2d.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart2d/piecewise_2D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart2d/bezier_2D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart2d/geometry.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart2d/solvelong.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart2d/piecestokes_2D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart2d/project.c)
rh_unit_test(test_formal_3d GEOMETRY cart3d TOPOLOGY THREE_D_PLANE
  test_formal_3d.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/piecewise_3D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/bezier_3D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/geometry.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/solvelong.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/piecestokes_3D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/project.c)
rh_unit_test(test_longchar_3d GEOMETRY cart3d TOPOLOGY THREE_D_PLANE
  test_longchar_3d.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/piecewise_3D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/bezier_3D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/geometry.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/solvelong.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/piecestokes_3D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/cart3d/project.c)
rh_unit_test(test_formal_stokes test_formal_stokes.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/piecewise_1D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
rh_unit_test(test_zeeman test_zeeman.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/bezier_1D.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/piecewise_1D.c)
rh_unit_test(test_ltepops test_ltepops.c)
rh_unit_test(test_ohchbf_range test_ohchbf_range.c)
rh_unit_test(test_molecule_range test_molecule_range.c)
rh_unit_test(test_error test_error.c)
rh_unit_test(test_status test_status.c)
rh_unit_test(test_readatom test_readatom.c
  ${PROJECT_SOURCE_DIR}/src/geometry/plane1d/project.c)
