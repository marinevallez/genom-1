#CMake usable only with specified version
cmake_minimum_required(VERSION 2.6)


#Project name : Genom1
project(Genom1 CXX C)

#’ON’ activates the tests; we can therefore use them.
option(test "Build tests." ON)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -std=c++11")
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_CXX_FLAGS "-W -Wall -pedantic -std=c++11")

include_directories(src)
link_directories("/usr/local/lib")


#Here we indicate the executables of interest
add_executable(Main src/main.cpp src/sequence.cpp src/matrixprotein.cpp src/Logo.cpp src/Letters.cpp src/Axes.cpp src/utilities.cpp)



set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH})
find_package(SFML 2 REQUIRED system window graphics network audio)
include_directories(${SFML_INCLUDE_DIR})
target_link_libraries(Main ${SFML_LIBRARIES})


#Generate Doxygen documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
	add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile 
				WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
				COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)


if (test)
#enable_testing()
#find_package(GTest)
#include_directories(${GTEST_INCLUDE_DIRS})
#add_executable(Genom1Test test/Genom1Test.cpp src/sequence.cpp src/#matrixprotein.cpp src/Logo.cpp src/Letters.cpp src/Axes.cpp)
#target_link_libraries(Genom1Test ${GTEST_BOTH_LIBRARIES} pthread)
#add_test(Genom1 Genom1_test)
endif(test)
