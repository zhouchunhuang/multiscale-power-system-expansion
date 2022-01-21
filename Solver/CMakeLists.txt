cmake_minimum_required(VERSION 3.15)
project(Model)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)

include_directories(inc)
include_directories(src)
include_directories(/usr/local/include)
include_directories(/Library/gurobi911/mac64/include)

# find GUROBI library
find_library(gurobi_static NAMES gurobi_g++4.2
        PATHS "/Library/gurobi911/mac64/lib")
find_library(gurobi_dynamic NAMES gurobi91
        PATHS "/Library/gurobi911/mac64/lib")
if(NOT gurobi_static)
    message(FATAL_ERROR "Failed to find Gurobi lib")
endif()
message(${gurobi_static})
message(${gurobi_dynamic})

# MPI library
find_library(mpi_lib NAMES mpi PATHS "/usr/local/lib")

add_executable(${PROJECT_NAME}
        src/main.cpp
        src/parameters.cpp
        src/Benders.cpp
        src/ColGen.cpp
        src/direct.cpp
        src/model.cpp
        src/NestColGen.cpp
        src/read.cpp
        src/test.cpp
        )
target_link_libraries(${PROJECT_NAME} ${gurobi_static} ${gurobi_dynamic} ${mpi_lib})

