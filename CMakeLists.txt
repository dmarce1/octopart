# Copyright (c) 2015 Thomas Heller
# Copyright (c) 2015 Dominic Marcello
# Copyright (c) 2018-2019 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.12.4)
cmake_policy(VERSION 3.12.4)

# Libraries linked via full path no longer produce linker search paths.
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

project(octopart CXX)

# Search path for CMake modules to be loaded by include() and find_package()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(HPX REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
find_package(Silo REQUIRED)
find_package(Boost REQUIRED)

set(common_source_files
	src/particle.cpp
	src/rand.cpp
)

set(octopart_source_files
	src/main.cpp
	src/math.cpp
	src/particle.cpp
	src/range.cpp
	src/state.cpp
	src/tree.cpp
	src/tree_stats.cpp
)

set(check2silo_source_files
    src/silo.cpp
    src/delaunay.cpp
)

set(common_header_files
	octopart/dim.hpp
	octopart/particle.hpp
	octopart/real.hpp
	octopart/state.hpp
	octopart/vect.hpp
)

set(octopart_header_files
	octopart/range.hpp
	octopart/rand.hpp
	octopart/tree.hpp
	octopart/tree_stats.hpp
)

set(check2silo_header_files
	octopart/delaunay.hpp	
)


add_hpx_executable(
  octopart
  DEPENDENCIES
    Silo::silo Boost::boost
  SOURCES
    ${octopart_source_files}
    ${common_source_files}
  HEADERS
    ${octopart_header_files}
    ${common_header_files}
)

add_hpx_executable(
  check2silo
  DEPENDENCIES
    Silo::silo Boost::boost
  SOURCES
    ${check2silo_source_files}
    ${common_source_files}
  HEADERS
    ${check2silo_header_files}
    ${common_header_files}
)

target_include_directories(octopart PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(check2silo PUBLIC ${PROJECT_SOURCE_DIR})


