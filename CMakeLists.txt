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

# Copy macros to build dir
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS mm_sim DESTINATION bin)
