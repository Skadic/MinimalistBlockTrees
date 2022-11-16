cmake_minimum_required(VERSION 3.0)
project(MinimalistBlockTrees)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O9")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")


add_subdirectory(external)

file(GLOB_RECURSE SOURCES "src/*.cpp")
add_library(blocktree STATIC ${SOURCES})
target_include_directories(blocktree PUBLIC includes)
target_include_directories(blocktree PRIVATE ${sdsl_SOURCE_DIR}/include)
target_link_libraries(blocktree INTERFACE sdsl)

# Add only if standalone
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    add_subdirectory(tests)
    add_subdirectory(executables)
endif()
