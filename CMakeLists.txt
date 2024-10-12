cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fno-inline-small-functions -DONIAK_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-g -O3")
set(CMAKE_CXX_STANDARD 23)
set(INCLUDE_TESTS ON)

project(doro)
find_package(libONIAK REQUIRED)

add_subdirectory(external/bch_codec)

add_executable(two_party "src/two_party.cpp")
include_directories("external")
include_directories("src")
target_link_libraries(two_party PRIVATE ONIAK::libONIAK)
target_link_libraries(two_party PRIVATE bch_codec)

if (INCLUDE_TESTS)
	add_subdirectory(test)
endif()

# add_executable(performance "src/performance.cpp")
# target_link_libraries(performance PRIVATE ONIAK::libONIAK)