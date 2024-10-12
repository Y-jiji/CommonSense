set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")

add_executable(rans_test "rans_test.cpp")
target_link_libraries(rans_test PRIVATE ONIAK::libONIAK)

add_executable(bch_test "bch_test.cpp")
target_link_libraries(bch_test PRIVATE ONIAK::libONIAK)
target_link_libraries(bch_test PRIVATE bch_codec)