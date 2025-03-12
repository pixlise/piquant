cmake_minimum_required(VERSION 3.10)

set(GIT_BRANCH "$ENV{CI_COMMIT_BRANCH}")

if ("${GIT_BRANCH}" STREQUAL "")
    set(GIT_BRANCH "$ENV{GIT_BRANCH_GIT}")
endif()

string(STRIP "${GIT_BRANCH}" GIT_BRANCH)

IF (WIN32)
    # Assume local build on windows, so we don't care about the branch
ELSE()
  # set stuff for other systems
  if ("${GIT_BRANCH}" STREQUAL "")
    message(FATAL_ERROR "git branch being built could not be determined")
  endif()
ENDIF()

message("CMake running for branch: ${GIT_BRANCH}")
message("CMake env var: BUILD_VERSION=$ENV{BUILD_VERSION}")

#if (NOT DEFINED "$ENV{BUILD_VERSION}")
#    message("BUILD_VERSION env var not defined, defining...")
#    set(ENV{BUILD_VERSION} "0.0.0")
#endif()

if ("$ENV{BUILD_VERSION}" STREQUAL "")
    message("BUILD_VERSION env var is empty, defining...")
    set(ENV{BUILD_VERSION} "0.0.0")
endif()

message("CMake project version used: $ENV{BUILD_VERSION}")

# set the project name and version
project(Piquant VERSION "$ENV{BUILD_VERSION}")

# set up template for version header that gets generated at build time
configure_file(version.h.in version.h)

if (DEFINED ENV{PIQUANT_EXCLUDE_PROTOBUF})
    # list all files needed to be compiled
    file(GLOB SOURCES
        src/*.h
        src/*.cpp
        src/data-formats/*.h
        src/data-formats/*.cc
    )
else()
    file(GLOB SOURCES
        src/*.h
        src/*.cpp
    )
endif()

# define executable to compile (using SOURCES defined above)
add_executable(Piquant ${SOURCES})

set_target_properties(Piquant PROPERTIES
                      ENABLE_EXPORTS 1
                      )

# include the directory that will contain version.h
target_include_directories(Piquant PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )

                           
if (DEFINED ENV{PIQUANT_EXCLUDE_PROTOBUF})
    message("PIQUANT_EXCLUDE_PROTOBUF not defined, so protobuf will be compiled in.")

    # Force static linking for protobuf, default was dynamic...
    find_library(PROTOBUF_LIB NAMES libprotobuf.a)
    #target_link_libraries(Piquant pthread protobuf-lite)
    target_link_libraries(Piquant pthread ${PROTOBUF_LIB})
else()
    message("PIQUANT_EXCLUDE_PROTOBUF is defined, so protobuf will be excluded.")
    
    add_compile_definitions(EXCLUDE_PIXLISE_READER)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set compiler flags to be pretty strict
IF (WIN32)
    # We set warning exclusions as command line parameters when running msbuild
ELSE()
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare")
ENDIF()
