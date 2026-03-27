cmake_minimum_required(VERSION 3.26)

project (THAMES_HYDRATION)

set (THAMES_HYDRATION_VERSION 5.2)
set (THAMES_HYDRATION_VERSION_MAJOR 5)
set (THAMES_HYDRATION_VERSION_MINOR 2)

# configure a header file to pass some CMake settings to source code

option(BUILD_DOC "Build API documentation" ON)
# look for doxygen capability
find_package (Doxygen)
if (DOXYGEN_FOUND)
  # set input and output files
  set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/src/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_SOURCE_DIR}/src/Doxyfile)
  # request to configure the file
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT}$ @ONLY)
  message("Doxygen build started")

  add_custom_target (doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

find_program (
  DOT
  NAMES dot
  HINTS "/opt/homebrew/bin/dot /opt/local/bin /usr/local/bin /usr/bin"
)

if (DOT)
  message("DOT found... yes")
  set (HAVE_DOT YES)
else()
  message("DOT found... no")
  set (HAVE_DOT NO)
  set (NO_DOT_WARNING "@warning Install graphviz/dot (http://www.graphviz.org)
  to create graphs and diagrams in this documentation")
endif()

find_program (
  PDFLATEX
  NAMES pdflatex
  HINTS "/Library/TeX/texbin /usr/bin /usr/local/bin"
)

if (PDFLATEX)
  message("pdflatex found... yes")
else()
  message("pdflatex found... no")
endif()

# Find PNG library for native PNG writing (removes ImageMagick dependency)
find_package(PNG REQUIRED)
if(PNG_FOUND)
    message("Found PNG: ${PNG_LIBRARIES}")
    include_directories(${PNG_INCLUDE_DIRS})
endif()

find_library (
    GEMS3K_LIB
    NAMES GEMS3K-static libGEMS3K-static libGEMS3K-static.a GEMS3K-static.a # what to look for
    HINTS "${CMAKE_SOURCE_DIR}/src/Resources/lib" # where to look
)

# check if we found the GEMS3K library
message(STATUS "GEMS3K_LIB: [${GEMS3K_LIB}]")

if (GEMS3K_LIB)
    message ("Found GEMSK3K_LIB: [${GEMS3K_LIB}]")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${GEMS3K_LIB})
    message ("EXTRA_LIBS: " ${EXTRA_LIBS})
else ()
    message (FATAL_ERROR "Did not find lib GEMS3K-static")
endif ()

find_library (
    MATH_LIB
    NAMES m libm # what to look for
    HINTS "/usr/lib" # where to look
)

message(STATUS "MATH_LIB: [${MATH_LIB}]")

if (MATH_LIB)
    message ("Found MATH_LIB: [${MATH_LIB}]")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${MATH_LIB})
else ()
    message (FATAL_ERROR "Did not find MATH lib")
endif ()

set (THAMES_HYDRATION_INCLUDE_PATH "-I${CMAKE_SOURCE_DIR}/src -I/usr/local/include")

message ("c++ compiler is : ${CMAKE_CXX_COMPILER_ID}")

# Set the list of compiler flags for GNU compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  set (CMAKE_CXX_FLAGS "-O2 -std=c++17 -DIPMGEMPLUGIN -Wall -Wno-misleading-indentation -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated ${THAMES_HYDRATION_INCLUDE_PATH}")
  #set (CMAKE_CXX_FLAGS "-O2 -std=c++17 -DIPMGEMPLUGIN ${THAMES_HYDRATION_INCLUDE_PATH}")
endif()

# Set the list of compiler flags for Clang compiler
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  add_compile_options(-O2 -std=c++17 -DIPMGEMPLUGIN -Wall -Wno-ignored-attributes -Wno-pedantic -Wno-variadic-macros -Wno-deprecated ${THAMES_HYDRATION_INCLUDE_PATH})
endif()

message ("Top level CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message ("Linked libraries are: " ${EXTRA_LIBS})

add_subdirectory (${CMAKE_SOURCE_DIR}/src/thameslib)

file (GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cc")

#set (EXTRA_LIBS ${EXTRA_LIBS} z)
#set (EXTRA_LIBS ${EXTRA_LIBS} m)

add_executable (thames ${SOURCES})

target_link_libraries (thames thameslib ${EXTRA_LIBS} PNG::PNG)

# macOS 26.4+ requires explicit ad-hoc signing for binaries to run
if(APPLE)
  add_custom_command(TARGET thames POST_BUILD
    COMMAND codesign --force --sign - $<TARGET_FILE:thames>
    COMMENT "Ad-hoc signing thames binary for macOS"
  )
endif()

install (TARGETS thames DESTINATION ${CMAKE_SOURCE_DIR}/bin)
if (EXISTS ${CMAKE_SOURCE_DIR}/build/html)
  install (DIRECTORY ${CMAKE_SOURCE_DIR}/build/html DESTINATION ${CMAKE_SOURCE_DIR}/doc)
endif ()
