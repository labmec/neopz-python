# This project uses pybind11 to create python binding for the NeoPZ library.
# Author: Gustavo Batistela
# Date: 13/09/2019

cmake_minimum_required (VERSION 3.11.0)

set (CMAKE_VERBOSE_MAKEFILE ON)

project (neopz-python)

# NeoPZ
IF (WIN32)
  # Prevents timespec redefinition problem with Visual Studio 2015
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
  find_package(PZ PATHS "C:/Program Files/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")
  find_package(PZ REQUIRED HINTS /home/natalia/CLionProjects/pzlibrary/neopzclean-install/CMake) 
ENDIF()

include_directories(${PZ_INCLUDE_DIRS})

# pybind11
find_package(pybind11)

# Project related sources
set(SOURCE_DIR "src")
pybind11_add_module(neopz ${SOURCES} "${SOURCE_DIR}/neopz_bindings.cpp")
target_link_libraries(neopz PRIVATE ${PZ_LIBRARIES})
