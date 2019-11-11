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
  find_package(PZ REQUIRED HINTS /Users/manouchehr/Documents/GitHub/neopz_python_install/CMake) 
ENDIF()

include_directories(${PZ_INCLUDE_DIRS})

# pybind11
find_package(pybind11)

# Project related sources
set(SOURCE_DIR "src")
pybind11_add_module(neopz ${SOURCES} "${SOURCE_DIR}/neopz_bindings.cpp")

IF(${PZ_USING_MKL})
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing -Wl,--no-as-needed")
  #find_library(MKL_LIB_DEF NAMES libmkl_def.so libmkl_def.dylib PATHS ${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
  #target_link_libraries(neopz PRIVATE ${PZ_LIBRARIES} ${MKL_LIB_DEF})
  target_link_libraries(neopz PRIVATE ${PZ_LIBRARIES})

ELSE()
  target_link_libraries(neopz PRIVATE ${PZ_LIBRARIES})
ENDIF()

