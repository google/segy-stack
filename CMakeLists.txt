# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required (VERSION 3.9)
cmake_policy(SET CMP0048 NEW)

project (segystack VERSION 0.1)

set(CMAKE_MACOSX_RPATH FALSE)
set(CMAKE_CXX_STANDARD 11)

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra" )
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g)
elseif(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    add_compile_options(-g -O0 --coverage)
    add_link_options(--coverage)
else()
    add_compile_options(-O3)
endif()

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

find_package(Protobuf REQUIRED)

file( GLOB_RECURSE segystack_HEADERS src/*.h )
file( GLOB_RECURSE segystack_SOURCES src/*.cc )
file( GLOB_RECURSE segystack_PROTOS src/*.proto )

include_directories(${Protobuf_INCLUDE_DIRS})
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${segystack_PROTOS})

add_library(${CMAKE_PROJECT_NAME} STATIC ${segystack_SOURCES} ${PROTO_SRCS} ${PROTO_HDRS})

add_subdirectory(util)
add_subdirectory(tests)
add_subdirectory(third_party/abseil-cpp)
add_subdirectory(third_party/glog)
add_subdirectory(third_party/pybind11)
add_subdirectory(python)

add_library(${CMAKE_PROJECT_NAME}_shared SHARED ${segystack_SOURCES} ${PROTO_SRCS} ${PROTO_HDRS})
set_target_properties( ${CMAKE_PROJECT_NAME}_shared PROPERTIES OUTPUT_NAME segystack )
set_target_properties( ${CMAKE_PROJECT_NAME}_shared PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION} )
target_link_libraries(${CMAKE_PROJECT_NAME}_shared ${Protobuf_LIBRARIES} glog::glog)

install( TARGETS ${CMAKE_PROJECT_NAME} ${CMAKE_PROJECT_NAME}_shared
         RUNTIME DESTINATION bin
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib
         PUBLIC_HEADER DESTINATION include/${libpath} )
