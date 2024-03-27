# This is the CMake script for compiling a CGAL application.

project( 2Dperiodic_ )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

# CGAL and its components
SET(CMAKE_BUILD_TYPE "Debug")
#ET(CMAKE_BUILD_TYPE "Release")
SET(CGAL_DIR $ENV{OWNPATH}/CGAL/cgal-4.4-${CMAKE_BUILD_TYPE}/lib/CGAL)
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++11 -g3")
  #ET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs -Wno-unused-parameter")
  #ET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-access-control")
  SET(CMAKE_VERBOSE_MAKEFILE false)
  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "available_volumes.cpp" )
  create_single_source_cgal_program( "available_volumes_after_takeout.cpp" )
  create_single_source_cgal_program( "unit_tests.cpp" )

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

