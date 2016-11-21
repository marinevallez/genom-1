cmake_minimum_required(VERSION 2.6)
project(Genom1 CXX C)
set(CMAKE_CXX_FLAGS "-W -Wall -pedantic -std=c++11")
add_executable(Main src/main.cpp src/sequence.cpp src/matrixprotein.cpp src/Logo.cpp src/Axes.cpp src/Letters.cpp)

set(CMAKE_MODULE_PATH "./cmake_modules" ${CMAKE_MODULE_PATH})
find_package(SFML 2 REQUIRED system window graphics network audio)
include_directories(${SFML_INCLUDE_DIR})
target_link_libraries(Main ${SFML_LIBRARIES})

find_package(Doxygen) 
add_custom_target(doc
	${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile WORKING_DIRECTORY $			{CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM)
