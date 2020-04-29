
file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cc)
set(SOURCES ${TEST_SOURCES})

set(TEST_BINARY ${CMAKE_PROJECT_NAME}_tests)

add_executable(${TEST_BINARY} ${TEST_SOURCES})
add_test(NAME ${TEST_BINARY} COMMAND ${TEST_BINARY})

target_link_libraries(${TEST_BINARY} PUBLIC ${CMAKE_PROJECT_NAME} gtest_main ${Protobuf_LIBRARIES} glog::glog)
