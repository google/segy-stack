
file(GLOB_RECURSE UTIL_SOURCES LIST_DIRECTORIES false *.h *.cc)

set(UTIL_BINARY ${CMAKE_PROJECT_NAME}_convert)

add_executable(${UTIL_BINARY} ${UTIL_SOURCES})

target_link_libraries(${UTIL_BINARY} PUBLIC ${CMAKE_PROJECT_NAME} absl::flags absl::flags_parse ${Protobuf_LIBRARIES} glog::glog )
