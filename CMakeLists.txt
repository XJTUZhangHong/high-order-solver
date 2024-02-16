cmake_minimum_required(VERSION 3.28)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})


SET(CMAKE_BUILD_TYPE "Release")
aux_source_directory(src SRC_SUB)
aux_source_directory(. SRC_CUR)
add_executable(solver ${SRC_SUB} ${SRC_CUR})
include_directories(include)
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3")
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

