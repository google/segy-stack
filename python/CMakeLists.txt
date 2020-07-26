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

find_package(PythonInterp)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

file(GLOB PYMOD_HEADERS "*.h")
file(GLOB PYMOD_SOURCES "*.cc")

pybind11_add_module( ${CMAKE_PROJECT_NAME}_py ${PYMOD_HEADERS} ${PYMOD_SOURCES} )

include_directories("${PROJECT_SOURCE_DIR}/third_party")

target_link_libraries( ${CMAKE_PROJECT_NAME}_py
                       PUBLIC ${CMAKE_PROJECT_NAME} ${Protobuf_LIBRARIES} glog::glog ${PYTHON_LIBRARIES})

set_target_properties( ${CMAKE_PROJECT_NAME}_py PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME}
                       PREFIX "${PYTHON_MODULE_PREFIX}"
                       SUFFIX ".so")

execute_process(COMMAND ${PYTHON_EXECUTABLE} -m site --user-site
                OUTPUT_VARIABLE PYTHON_USERSITE
                OUTPUT_STRIP_TRAILING_WHITESPACE)

install(FILES $<TARGET_FILE:${CMAKE_PROJECT_NAME}_py> DESTINATION ${PYTHON_USERSITE})
