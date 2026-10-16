# RH: multi-level, multi-dimensional NLTE radiative transfer.
#
#   cmake -B build                 configure (default build type RelWithDebInfo)
#   cmake --build build -j         build librh, the four programs and the tools
#   ctest --test-dir build         quick check that the programs start
#   cmake --install build --prefix ~/rh
#
# Options (-D<option>=<value>):
#   RH_NATIVE=ON        tune for the build machine (-march=native or
#                       -mcpu=native); the binaries may not run elsewhere
#   RH_SANITIZE=<list>  build with sanitizers, e.g. address or
#                       address,undefined
#   RH_WARNINGS=ON      compile with -Wall -Wextra
#   RH_BUILD_TOOLS=OFF  skip the atomic-data tools in tools/

cmake_minimum_required(VERSION 3.20)
project(RH VERSION 3.0 LANGUAGES C)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Build type: Debug, Release, RelWithDebInfo or MinSizeRel" FORCE)
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(RH_NATIVE "Tune for the build machine's CPU" OFF)
option(RH_WARNINGS "Compile with -Wall -Wextra" OFF)
option(RH_BUILD_TOOLS "Build the atomic-data tools in tools/" ON)
set(RH_SANITIZE "" CACHE STRING
    "Comma-separated sanitizers, e.g. address or address,undefined")

include(CheckCCompilerFlag)
include(CheckIncludeFile)
include(CheckSymbolExists)

# --- Compiler settings shared by every target ----------------------------

add_library(rh_options INTERFACE)

# RH relies on IEEE semantics (exact comparisons in the wavelength grid,
# FP exception traps). Intel compilers default to value-changing
# optimizations, so ask for precise floating point there.
if(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
  target_compile_options(rh_options INTERFACE -fp-model=precise)
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  target_compile_options(rh_options INTERFACE "SHELL:-fp-model precise")
endif()

if(RH_NATIVE)
  check_c_compiler_flag(-march=native RH_HAVE_MARCH_NATIVE)
  check_c_compiler_flag(-mcpu=native RH_HAVE_MCPU_NATIVE)
  if(RH_HAVE_MARCH_NATIVE)
    target_compile_options(rh_options INTERFACE -march=native)
  elseif(RH_HAVE_MCPU_NATIVE)
    target_compile_options(rh_options INTERFACE -mcpu=native)
  else()
    message(WARNING "RH_NATIVE: the compiler accepts neither -march=native "
                    "nor -mcpu=native; building generic code")
  endif()
endif()

if(RH_WARNINGS)
  target_compile_options(rh_options INTERFACE -Wall -Wextra)
endif()

if(RH_SANITIZE)
  set(san -fsanitize=${RH_SANITIZE} -fno-omit-frame-pointer)
  target_compile_options(rh_options INTERFACE ${san})
  target_link_options(rh_options INTERFACE ${san})
endif()

# --- Dependencies ----------------------------------------------------------

find_package(Threads REQUIRED)
find_library(RH_MATH_LIBRARY m)

# XDR: newer glibc no longer provides Sun RPC, so Linux systems use
# libtirpc; macOS and some other systems still have XDR in the C library.
add_library(rh_xdr INTERFACE)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(TIRPC QUIET IMPORTED_TARGET libtirpc)
endif()
if(TIRPC_FOUND)
  target_link_libraries(rh_xdr INTERFACE PkgConfig::TIRPC)
  message(STATUS "XDR: libtirpc ${TIRPC_VERSION}")
else()
  find_path(TIRPC_INCLUDE_DIR tirpc/rpc/xdr.h)
  if(TIRPC_INCLUDE_DIR)
    set(TIRPC_INCLUDE_DIR ${TIRPC_INCLUDE_DIR}/tirpc)
  endif()
  find_library(TIRPC_LIBRARY tirpc)
  if(TIRPC_INCLUDE_DIR AND TIRPC_LIBRARY)
    target_include_directories(rh_xdr INTERFACE ${TIRPC_INCLUDE_DIR})
    target_link_libraries(rh_xdr INTERFACE ${TIRPC_LIBRARY})
    message(STATUS "XDR: libtirpc (${TIRPC_LIBRARY})")
  else()
    check_symbol_exists(xdr_double "rpc/types.h;rpc/xdr.h" RH_HAVE_LIBC_XDR)
    if(NOT RH_HAVE_LIBC_XDR)
      message(FATAL_ERROR "No XDR library found. Install libtirpc and its "
              "headers (libtirpc-devel on RHEL/Fedora/SUSE, libtirpc-dev "
              "on Debian/Ubuntu).")
    endif()
    message(STATUS "XDR: C library")
  endif()
endif()

# Floating-point exception traps and backtraces (fpehandler.c)
set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
if(RH_MATH_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${RH_MATH_LIBRARY})
endif()
check_symbol_exists(feenableexcept fenv.h RH_HAVE_FEENABLEEXCEPT)
unset(CMAKE_REQUIRED_DEFINITIONS)
unset(CMAKE_REQUIRED_LIBRARIES)
check_include_file(execinfo.h RH_HAVE_EXECINFO_H)

# --- The library -----------------------------------------------------------

add_library(rh STATIC
  abundance.c accelerate.c background.c backgropac_xdr.c barklem.c
  bezier_aux.c broad.c brs_xdr.c chemequil.c cocollisions.c collision.c
  cubeconvol.c duplicate.c error.c expint.c expspline.c fillgamma.c
  fixedrate.c fpehandler.c gammafunc.c gaussleg.c getcpu.c getlambda.c
  getline.c giigen.c h2collisions.c humlicek.c hunt.c hydrogen.c
  initial_xdr.c initscatter.c iterate.c kurucz.c linear.c ltepops.c
  ludcmp.c matrix.c maxchange.c metal.c molzeeman.c nemetals.c ohchbf.c
  opacity.c options.c order.c parallel.c parse.c paschen.c planck.c
  pops_xdr.c profile.c radrate_xdr.c rayleigh.c readatom.c readb_xdr.c
  readinput.c readj.c readmolecule.c readvalue.c redistribute.c scatter.c
  solvene.c sortlambda.c spline.c statequil.c stokesopac.c stopreq.c
  thomson.c vacuumtoair.c voigt.c w3.c wigner.c writeatmos_xdr.c
  writeatom_xdr.c writecoll_xdr.c writedamp_xdr.c writeinput_xdr.c
  writemetal_xdr.c writemolec_xdr.c writeopac_xdr.c writespect_xdr.c
  zeeman.c)
target_include_directories(rh PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(rh PUBLIC rh_options rh_xdr Threads::Threads)
if(RH_MATH_LIBRARY)
  target_link_libraries(rh PUBLIC ${RH_MATH_LIBRARY})
endif()
if(RH_HAVE_FEENABLEEXCEPT)
  set_property(SOURCE fpehandler.c APPEND PROPERTY COMPILE_DEFINITIONS
               HAVE_FEENABLEEXCEPT)
endif()
if(RH_HAVE_EXECINFO_H)
  set_property(SOURCE fpehandler.c APPEND PROPERTY COMPILE_DEFINITIONS
               HAVE_EXECINFO_H)
endif()
# Complex arithmetic in the Voigt routines without the NaN/infinity
# recovery of C99 Annex G, which the arguments there never need: faster
# multiplication. (Division uses cdiv() from cdiv.h, so results do not
# depend on this flag; Apple clang, for one, does not have it.)
check_c_compiler_flag(-fcx-limited-range RH_HAVE_CX_LIMITED_RANGE)
if(RH_HAVE_CX_LIMITED_RANGE)
  set_property(SOURCE voigt.c humlicek.c APPEND PROPERTY COMPILE_OPTIONS
               -fcx-limited-range)
endif()

# --- Programs --------------------------------------------------------------

# rh_program(TARGET NAME SOURCES...): an executable called NAME in the
# build directory of the current geometry (build/rhf1d/solveray, ...). It
# is installed as INSTALL_NAME (default NAME), since several geometries
# have programs of the same name.
function(rh_program target name)
  cmake_parse_arguments(P "" "INSTALL_NAME" "" ${ARGN})
  add_executable(${target} ${P_UNPARSED_ARGUMENTS})
  set_target_properties(${target} PROPERTIES OUTPUT_NAME ${name})
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${target} PRIVATE rh)
  if(NOT P_INSTALL_NAME)
    set(P_INSTALL_NAME ${name})
  endif()
  install(PROGRAMS $<TARGET_FILE:${target}> TYPE BIN RENAME ${P_INSTALL_NAME})
endfunction()

add_subdirectory(rhf1d)
add_subdirectory(rhsc2d)
add_subdirectory(rhsc3d)
add_subdirectory(rhsphere)
if(RH_BUILD_TOOLS)
  add_subdirectory(tools)
endif()

# --- Smoke tests -----------------------------------------------------------
# Each program, run where it finds no input, must fail cleanly with RH's
# own error message rather than crash. The physics regression suite is
# tests/run_tests.sh.

enable_testing()
foreach(prog rhf1d rhsc2d rhsc3d rhsphere)
  add_test(NAME start_${prog} COMMAND ${prog}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(start_${prog} PROPERTIES
    PASS_REGULAR_EXPRESSION "TERMINATING_ERROR"
    FAIL_REGULAR_EXPRESSION "Sanitizer|runtime error|Segmentation|FATAL ERROR")
endforeach()

message(STATUS "RH ${PROJECT_VERSION}: build type ${CMAKE_BUILD_TYPE}, "
               "native=${RH_NATIVE}, sanitize='${RH_SANITIZE}', "
               "warnings=${RH_WARNINGS}")
