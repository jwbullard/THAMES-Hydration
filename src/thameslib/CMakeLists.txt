cmake_minimum_required(VERSION 3.26)

project (ThamesLibrary_project)

file (GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/thameslib/*.cc")

message ("THAMES lib directory, CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

add_library (thameslib STATIC ${SOURCES})

# target_link_libraries(thameslib ${EXTRA_LIBS})
target_link_libraries(thameslib)

install (TARGETS thameslib DESTINATION ${CMAKE_SOURCE_DIR}/lib)
