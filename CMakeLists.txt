cmake_minimum_required(VERSION 2.6)
project(Genom1 CXX C)
set(CMAKE_CXX_FLAGS "-W -Wall -pedantic -std=c++11")
add_executable(Main src/main.cc)
target_link_libraries( )


find_package(Doxygen) 
add_custom_target(doc