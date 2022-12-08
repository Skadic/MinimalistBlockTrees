cmake_minimum_required(VERSION 3.0)
project(MinimalistBlockTrees)

set(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O9")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp")

add_compile_options (-fdiagnostics-color=always)


add_subdirectory(external)

file(GLOB_RECURSE SOURCES "src/*.cpp")
add_library(blocktree STATIC ${SOURCES})
target_include_directories(blocktree PUBLIC includes)
target_include_directories(blocktree SYSTEM PRIVATE ${sdsl_SOURCE_DIR}/include)
set_target_properties(sdsl PROPERTIES COMPILE_FLAGS "-w")
target_link_libraries(blocktree INTERFACE sdsl)

# Add only if standalone
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    add_subdirectory(tests)
    add_subdirectory(executables)
endif()
