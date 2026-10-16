cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(MicromegasSim VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Geant4 -- all physics modules
find_package(Geant4 REQUIRED ui_all vis_all)

# Use Geant4 include dirs and definitions
include(${Geant4_USE_FILE})

# Find ROOT (optional but recommended for output)
find_package(ROOT QUIET)

# Include dirs
include_directories(${PROJECT_SOURCE_DIR}/include)

# ── Source revision stamped into the binary ──────────────────────────────
# Every output file records the commit it was produced by (RunMeta.cc ->
# the RunMeta tree), so a ROOT file separated from its job log is still
# identifiable. SIM_CAMPAIGN_PLAN §9: no un-manifested runs.
#
# This resolves at CONFIGURE time, so .git/HEAD and the index are added to
# CMAKE_CONFIGURE_DEPENDS -- committing or checking out re-runs cmake and
# restamps. A dirty tree is recorded as such rather than hidden; for a
# production run that flag is itself a finding.
find_package(Git QUIET)
set(P2_GIT_HASH "unknown")
set(P2_GIT_DIRTY 0)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE P2_GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE _p2_git_status
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(NOT "${_p2_git_status}" STREQUAL "")
        set(P2_GIT_DIRTY 1)
    endif()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
                 "${PROJECT_SOURCE_DIR}/.git/HEAD"
                 "${PROJECT_SOURCE_DIR}/.git/index")
endif()
message(STATUS "Source revision: ${P2_GIT_HASH} (dirty: ${P2_GIT_DIRTY})")

# Source files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)

# ── Garfield++ (optional): Heed for the drift-gas ionization, --heed ─────
# Must be built against the SAME ROOT and Geant4 as this program. On lxplus
# that is the Garfield++ of the LCG view used by setup_lxplus.sh (LCG_106:
# Geant4 11.2.1, ROOT 6.32.02, commit d3bbb29); its Heed was checked against
# the pinned standalone build (same spectra, docs/results/heed/pai).
set(P2_GARFIELD_DIR "$ENV{GARFIELD_INSTALL}" CACHE PATH "Garfield++ install prefix")
if(NOT P2_GARFIELD_DIR)
    set(P2_GARFIELD_DIR "/cvmfs/sft.cern.ch/lcg/releases/LCG_106/Garfield++/HEAD/x86_64-el9-gcc13-opt")
endif()
find_path(GARFIELD_INCLUDE_DIR Garfield/TrackHeed.hh HINTS ${P2_GARFIELD_DIR}/include NO_DEFAULT_PATH)
find_library(GARFIELD_LIBRARY Garfield HINTS ${P2_GARFIELD_DIR}/lib64 ${P2_GARFIELD_DIR}/lib NO_DEFAULT_PATH)
set(P2_WITH_GARFIELD OFF)
if(GARFIELD_INCLUDE_DIR AND GARFIELD_LIBRARY AND ROOT_FOUND)
    set(P2_WITH_GARFIELD ON)
    message(STATUS "Garfield++ found: ${GARFIELD_LIBRARY} (--heed enabled)")
else()
    list(FILTER SOURCES EXCLUDE REGEX "HeedGasModel\\.cc$")
    message(STATUS "Garfield++ not found under ${P2_GARFIELD_DIR}: --heed disabled")
endif()

# Build the executable
add_executable(mm_sim mm_sim.cc ${SOURCES})
target_link_libraries(mm_sim ${Geant4_LIBRARIES})
target_compile_definitions(mm_sim PRIVATE
    P2_GIT_HASH="${P2_GIT_HASH}"
    P2_GIT_DIRTY=${P2_GIT_DIRTY})

if(ROOT_FOUND)
    target_include_directories(mm_sim PRIVATE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(mm_sim ${ROOT_LIBRARIES})
    target_compile_definitions(mm_sim PRIVATE USE_ROOT)
    message(STATUS "ROOT found and enabled: ${ROOT_VERSION}")
else()
    message(STATUS "ROOT not found; CSV output only")
endif()

if(P2_WITH_GARFIELD)
    target_include_directories(mm_sim PRIVATE ${GARFIELD_INCLUDE_DIR})
    target_link_libraries(mm_sim ${GARFIELD_LIBRARY})
    get_filename_component(_garf_libdir ${GARFIELD_LIBRARY} DIRECTORY)
    set_target_properties(mm_sim PROPERTIES BUILD_RPATH "${_garf_libdir}" INSTALL_RPATH "${_garf_libdir}")
    target_compile_definitions(mm_sim PRIVATE P2_WITH_GARFIELD
        P2_HEED_DATABASE="${P2_GARFIELD_DIR}/share/Heed/database")
    # Older Garfield (LCG_106) has a global randomEngine, no Random::SetEngine.
    file(READ ${GARFIELD_INCLUDE_DIR}/Garfield/Random.hh _garf_random)
    string(FIND "${_garf_random}" "SetEngine" _garf_setengine)
    if(_garf_setengine EQUAL -1)
        target_compile_definitions(mm_sim PRIVATE GARFIELD_OLD_RANDOM)
    endif()
endif()

# Copy macros to build dir
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS mm_sim DESTINATION bin)
