cmake_minimum_required(VERSION 3.17)
project(apitest-proj)
add_executable(apitest apiTest2.cc)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
target_link_libraries(apitest PUBLIC ${MPI_CXX_LIBRARIES})


find_library(DFTFE_LIBRARY
    NAMES dftfem
    HINTS ${DFTFE_INSTALL_PATH}/lib
    NO_DEFAULT_PATH
  )
 
IF(DFTFE_LIBRARY STREQUAL "DFTFE_LIBRARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "-- DFTFE was not found in ${DFTFE_INSTALL_PATH}")
ENDIF()


target_include_directories(apitest PUBLIC "${DFTFE_INSTALL_PATH}/include")
target_link_libraries(apitest PUBLIC ${DFTFE_LIBRARY})
