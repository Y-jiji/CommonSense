cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 20)

project(doro)
find_package(libONIAK REQUIRED)

add_executable(doro_exp "src/doro_exp.cpp")
target_link_libraries(doro_exp PRIVATE ONIAK::libONIAK)