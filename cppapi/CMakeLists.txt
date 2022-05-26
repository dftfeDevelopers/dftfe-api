cmake_minimum_required(VERSION 3.17)
project(dftfeapitest-proj)

add_executable(dftfeapitest example1/apiTestSinglePeriodicSystemMPI.cpp)
#add_executable(dftfeapitest example2/apiTestTwoSystemsMPI.cpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
target_link_libraries(dftfeapitest PUBLIC ${MPI_CXX_LIBRARIES})


find_library(DFTFE_LIBRARY
    NAMES dftfem
    HINTS ${DFTFE_INSTALL_PATH}/lib
    NO_DEFAULT_PATH
  )
 
IF(DFTFE_LIBRARY STREQUAL "DFTFE_LIBRARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "-- DFTFE was not found in ${DFTFE_INSTALL_PATH}")
ENDIF()


target_include_directories(dftfeapitest PUBLIC "${DFTFE_INSTALL_PATH}/include")
target_link_libraries(dftfeapitest PUBLIC ${DFTFE_LIBRARY})
