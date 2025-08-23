cmake_minimum_required(VERSION 3.12)
project(CQOI VERSION 0.0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-g -O3")

include(GNUInstallDirs)
message(STATUS "CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")

find_package(ADIOS2 REQUIRED CONFIG
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/ADIOS2/adios2-install"
  PATH_SUFFIXES ${CMAKE_INSTALL_LIBDIR}/cmake/adios2
)

find_path(ZSTD_INCLUDE_DIR NAMES zstd.h
  HINTS
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/tools/zstd/lib"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/tools/zstd/lib/common"
    /usr/include /usr/local/include
)

find_library(ZSTD_LIBRARY NAMES zstd libzstd
  HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ2/install"
  PATH_SUFFIXES lib ${CMAKE_INSTALL_LIBDIR}
)

if(NOT ZSTD_INCLUDE_DIR OR NOT ZSTD_LIBRARY)
  message(FATAL_ERROR
    "Could not find zstd: header=${ZSTD_INCLUDE_DIR}, lib=${ZSTD_LIBRARY}\n"
    "Please ensure libzstd is installed or built in SZ2/SZ3.")
endif()

find_library(SZ3_LIBRARY NAMES SZ3 sz3
  HINTS
    "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/build"
  PATH_SUFFIXES lib src
)

set(SZ3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/SZ3/include")

if(NOT SZ3_LIBRARY)
  message(FATAL_ERROR "Could not find SZ3 library in external/SZ3/build/lib or build/src")
endif()

set (MGARDx_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/external/MGARDx/install/include")

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${SZ3_INCLUDE_DIR}"
  "${ZSTD_INCLUDE_DIR}"
)

target_link_libraries(${PROJECT_NAME} INTERFACE
  "${ZSTD_LIBRARY}"
  "${SZ3_LIBRARY}"
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

add_subdirectory(test)
add_subdirectory(parallel_src)
