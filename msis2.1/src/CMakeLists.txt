cmake_minimum_required(VERSION 3.15...3.21)

project(MSIS2
  LANGUAGES Fortran
  DESCRIPTION "Build and use MSIS 2.1 from Python"
  HOMEPAGE_URL https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.1/
)

set(FETCHCONTENT_UPDATES_DISCONNECTED_MSIS2 true)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

add_subdirectory(msis2.1)

add_executable(int_driver intdriver.f90)
target_link_libraries(int_driver PRIVATE msis2)

find_package (NetCDF COMPONENTS C Fortran REQUIRED)
include_directories(${NETCDF_INCLUDE_DIRS})
target_link_libraries (int_driver PRIVATE ${NETCDF_Fortran_LIBRARIES})

set_target_properties(int_driver PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# MSIS 2.1 exe needs this parm file
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/msis21.parm)
  file(COPY ${msis2proj_SOURCE_DIR}/msis21.parm DESTINATION ${PROJECT_SOURCE_DIR})
endif()
