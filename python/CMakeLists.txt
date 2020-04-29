file(GLOB PYMOD_HEADERS "*.h")
file(GLOB PYMOD_SOURCES "*.cc")

# set(PYBIND11_PYTHON_VERSION 2.7)

pybind11_add_module( ${CMAKE_PROJECT_NAME}_py ${PYMOD_HEADERS} ${PYMOD_SOURCES})

include_directories("${PROJECT_SOURCE_DIR}/third_party")

target_link_libraries( ${CMAKE_PROJECT_NAME}_py PUBLIC ${CMAKE_PROJECT_NAME}_shared )

set_target_properties( ${CMAKE_PROJECT_NAME}_py PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME}
                       PREFIX "${PYTHON_MODULE_PREFIX}"
                       SUFFIX ".so")
