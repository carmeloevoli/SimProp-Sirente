cmake_minimum_required(VERSION 3.1)
project(simprop CXX)
set(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(SIMPROP_EXTRA_SOURCES)
set(SIMPROP_EXTRA_INCLUDES)
set(SIMPROP_EXTRA_LIBRARIES)

set(CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0048 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -pedantic -ftree-vectorize -O3")
set(SIMPROP_INCLUDES)
set(SIMPROP_LIBRARIES)

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
    #set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# PLOG (provided)
include_directories(external/plog/include)

# NamedType (provided)
include_directories(external/NamedType/include)

# C++ Threads required
#find_package(Threads REQUIRED)
#list(APPEND SIMPROP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# GSL (required)
find_package(GSL REQUIRED)
list(APPEND SIMPROP_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
list(APPEND SIMPROP_EXTRA_LIBRARIES ${GSL_LIBRARIES})

# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
option(ENABLE_TESTING "Build tests and enable test target" OFF)
if(ENABLE_TESTING)
    set(GOOGLETEST_VERSION 1.10.0)
    include_directories(external/gtest/googletest/include)
    add_subdirectory(external/gtest/googletest)
    if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
    endif(APPLE)
endif(ENABLE_TESTING)

# Version info from Git
option(ENABLE_GIT "Embedding information about Simprop version from git" ON)
if(ENABLE_GIT)
    include(GetGitRevisionDescription)
    get_git_head_revision(GIT_REFSPEC GIT_SHA1)
    git_local_changes(GIT_HAS_LOCAL_CHANGES)
    set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
    configure_file("src/utils/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)
endif(ENABLE_GIT)
IF(NOT ENABLE_GIT OR (${GIT_SHA1} STREQUAL "GITDIR-NOTFOUND"))
    set(GIT_REFSPEC "")
    set(GIT_SHA1 "")
    set(GIT_DESC "${SIMPROP_RELEASE_VERSION}-no-git")
endif()
message(STATUS "Simprop version: ${GIT_DESC} ${GIT_SHA1} ${GIT_REFSPEC}")

# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
include_directories(include ${SIMPROP_EXTRA_INCLUDES})

add_library(simprop SHARED
    src/analyticalSolutions/beniamino.cpp
    src/crossSections/PhotoPionProductionXsec.cpp
    src/energyLosses/AdiabaticContinuousLosses.cpp
    src/energyLosses/BGG2002ContinuousLosses.cpp
    src/energyLosses/PairProductionLosses.cpp
    src/energyLosses/PhotoPionContinuousLosses.cpp
    src/interactions/PhotoPionProduction.cpp
    src/particleStacks/Builder.cpp
    src/particleStacks/SingleParticleBuilder.cpp
    src/particleStacks/SingleSourceBuilder.cpp
    src/particleStacks/SourceEvolutionBuilder.cpp
    src/photonFields/LookupTablePhotonField.cpp
    src/photonFields/CmbPhotonField.cpp
	src/utils/io.cpp
	src/utils/logging.cpp
    src/utils/numeric.cpp
    src/utils/timer.cpp
    src/common.cpp
    src/cosmology.cpp
   	#src/params.cpp
 	#src/particleStack.cpp
    src/pid.cpp
    "${git_revision_cpp}"
    )
    target_link_libraries(simprop ${SIMPROP_EXTRA_LIBRARIES})

# make library
#include_directories(include ${SIMPROP_INCLUDES})
#add_library(SIMPROP_LIB ${SIMPROP_SRCS} ${SIMPROP_INCLUDES})

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# copy data files in build
file(GLOB SIMPROP_DATA "data/*.txt")
file(COPY ${SIMPROP_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# make examples
add_executable(print_cosmology examples/printCosmology.cpp)
target_link_libraries (print_cosmology simprop ${SIMPROP_EXTRA_LIBRARIES})

add_executable(print_photons examples/printPhotonFields.cpp)
target_link_libraries (print_photons simprop ${SIMPROP_EXTRA_LIBRARIES})

add_executable(print_xsecs examples/printCrossSections.cpp)
target_link_libraries (print_xsecs simprop ${SIMPROP_EXTRA_LIBRARIES})

add_executable(print_losses examples/printLosses.cpp)
target_link_libraries (print_losses simprop ${SIMPROP_EXTRA_LIBRARIES})

add_executable(print_initial examples/printInitialState.cpp)
target_link_libraries (print_initial simprop ${SIMPROP_EXTRA_LIBRARIES})

add_executable(print_rates examples/printRates.cpp)
target_link_libraries (print_rates simprop ${SIMPROP_EXTRA_LIBRARIES})

add_executable(analytical_solution examples/analyticalSolution.cpp)
target_link_libraries (analytical_solution simprop ${SIMPROP_EXTRA_LIBRARIES})

add_executable(evolve examples/evolveContinuousLosses.cpp)
target_link_libraries (evolve simprop ${SIMPROP_EXTRA_LIBRARIES})

# add_executable(sandbox src/sandbox.cpp)
# target_link_libraries (sandbox SIMPROP_LIB ${SIMPROP_LIBRARIES})

# testing
if(ENABLE_TESTING)
    enable_testing()
    add_executable(testCommon test/testCommon.cpp)
    target_link_libraries(testCommon simprop gtest gtest_main ${SIMPROP_EXTRA_LIBRARIES})
    add_test(testCommon testCommon)

    add_executable(testCosmology test/testCosmology.cpp)
    target_link_libraries(testCosmology simprop gtest gtest_main ${SIMPROP_EXTRA_LIBRARIES})
    add_test(testCosmology testCosmology)

    add_executable(testRng test/testRng.cpp)
    target_link_libraries(testRng simprop gtest gtest_main ${SIMPROP_EXTRA_LIBRARIES})
    add_test(testRng testRng)

    add_executable(testPid test/testPid.cpp)
    target_link_libraries(testPid simprop gtest gtest_main ${SIMPROP_EXTRA_LIBRARIES})
    add_test(testPid testPid)
endif(ENABLE_TESTING)

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "SIMPROP_INCLUDES ...... = ${SIMPROP_EXTRA_INCLUDES}")
message (STATUS "SIMPROP_LIBRARIES ..... = ${SIMPROP_EXTRA_LIBRARIES}")
