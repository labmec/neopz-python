# This project uses pybind11 to create python binding for the NeoPZ library.
# Author: Gustavo Batistela
# Date: 13/09/2019

cmake_minimum_required (VERSION 3.11.0)

set (CMAKE_VERBOSE_MAKEFILE ON)

project (neopz-bindings)

# NeoPZ
IF (WIN32)
  # Prevents timespec redefinition problem with Visual Studio 2015
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
  find_package(PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")
  find_package(PZ REQUIRED)
ENDIF()

include_directories(${PZ_INCLUDE_DIRS})

# Locations to search for external libraries
set(SEARCH_DIRS ENV PATH
  /usr/local/opt/
  /usr/include /usr/local/include
  /usr/local/include/boost
  /usr/local/include/log4cxx
  /usr/lib /usr/lib64
  /usr/local/lib
  ./.
  ./..
  ./../libs
  ./../libs/fad
  ./../libs/metis
  ./../libs/metis/include
  ./../libs/metis/lib
  ./../externallibs
  ./../externallibs/pthread
  ./externallibs/lib
  ./externallibs/include
  ./../externallibs/lib
  ./../externallibs/include
)

# log4cxx
find_library(LOG4CXX_LIB_RELEASE
        NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a
        PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Release
        )
find_library(LOG4CXX_LIB_DEBUG
        NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a       # Change name here if debug lib has different name.
        PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Debug
        )
include_directories(${LOG4CXX_INCLUDE})

# PThread
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/include )
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} D:/Projetos/externallibs/lib  ${PROJECT_SOURCE_DIR}/../pthread/lib)

# Project related sources
set(SOURCE_DIR "src/neopz")
include_directories(${SOURCE_DIR})

find_package(pybind11)
pybind11_add_module(neopz ${SOURCES} "${SOURCE_DIR}/neopz_bindings.cpp")
target_link_libraries(neopz PRIVATE ${PZ_LIBRARIES})
