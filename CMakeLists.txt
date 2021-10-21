# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.7)
add_compile_options(-g)
# Project's name

project(electron_heat_capacity)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")

find_package(Python 3.9 REQUIRED)
find_package(PythonLibs REQUIRED) 
find_package(GSL REQUIRED)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
link_directories(/home/boyang/softwares/yaml/build)

include_directories(
	${PROJECT_SOURCE_DIR}/include
	${PYTHON_INCLUDE_DIRS}
	/home/boyang/softwares/yaml/yaml-cpp/include/
	/home/boyang/.local/lib/python3.9/site-packages/numpy/core/include
)
file(GLOB libsource "src/*")
file(GLOB libinclude "include/*")
add_library(plot SHARED ${libsource} ${libinclude})
target_link_libraries(plot ${PYTHON_LIBRARIES} ${GSL_LIBRARIES})

file( GLOB APP_SOURCES RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/util/* )
foreach( sourcefile ${APP_SOURCES} )
    string( REPLACE ".C" "_run" exename ${sourcefile} )
    string( REPLACE "util/" "" exename ${exename} )
    add_executable( ${exename} ${sourcefile})
	target_link_libraries(${exename} plot yaml-cpp)
endforeach( sourcefile ${APP_SOURCES} )
