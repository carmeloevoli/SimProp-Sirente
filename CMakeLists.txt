cmake_minimum_required(VERSION 3.1)
set(LANGS CXX)
project(simprop ${LANGS})
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION VERSION)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)
cmake_policy(SET CMP0048 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")
set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -pedantic -ftree-vectorize -O3 -DDEBUG")
set(SIMPROP_INCLUDES)
set(SIMPROP_LIBRARIES)

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    #set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
    message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)


# ----------------------------------------------------------------------------
# Fix Apple RPATH
# ----------------------------------------------------------------------------
set(CMAKE_MACOSX_RPATH 1)
option(USE_ABSOLUTE_RPATH "Add absolute rpath to all libraries and executables" ON)
if(APPLE OR USE_ABSOLUTE_RPATH)
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(ABSOLUTE_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
    if(NOT IS_ABSOLUTE ${ABSOLUTE_RPATH})
        set(ABSOLUTE_RPATH ${CMAKE_BINARY_DIR}/${ABSOLUTE_RPATH})
    endif()

    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${ABSOLUTE_RPATH}" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
        message(STATUS "Use absolute RPATH ${ABSOLUTE_RPATH}")
        set(CMAKE_INSTALL_RPATH "${ABSOLUTE_RPATH}")
    endif("${isSystemDir}" STREQUAL "-1")
endif()

# Fix for linker error on mac, Issue: #147
if(APPLE)
    add_definitions(-arch x86_64)
endif(APPLE)


# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# PLOG (provided)
include_directories(external/plog/include)

# C++ Threads required
find_package(Threads REQUIRED)
list(APPEND SIMPROP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# GSL required
find_package(GSL REQUIRED)
list(APPEND SIMPROP_INCLUDES ${GSL_INCLUDE_DIR})
list(APPEND SIMPROP_LIBRARIES ${GSL_LIBRARIES})

# Version info from Git
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_HAS_LOCAL_CHANGES)
set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
configure_file("src/utils/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------

include_directories("${CMAKE_SOURCE_DIR}")
include_directories(include)

set(SIMPROP_SRCS
    src/common.cpp
 	src/params.cpp
    src/pid.cpp
    src/simprop.cpp
    src/photonFields/CMB.cpp
	src/utils/logging.cpp
    src/utils/misc.cpp
    src/utils/timer.cpp
    "${git_revision_cpp}"
    )

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# make library
include_directories(include ${SIMPROP_INCLUDES})
add_library(SIMPROP_LIB ${SIMPROP_SRCS} ${SIMPROP_INCLUDES})

# make executable
add_executable(simprop src/main.cpp)
target_link_libraries (simprop SIMPROP_LIB ${SIMPROP_LIBRARIES})

add_executable(sandbox src/sandbox.cpp)
target_link_libraries (sandbox SIMPROP_LIB ${SIMPROP_LIBRARIES})

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "SIMPROP_INCLUDES ........ = ${SIMPROP_INCLUDES}")
message (STATUS "SIMPROP_LIBRARIES ....... = ${SIMPROP_LIBRARIES}")