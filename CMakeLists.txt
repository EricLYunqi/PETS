cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O3 -Wall -g2 -ggdb -fopenmp -mavx2")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -fopenmp -mavx2")
SET(CMAKE_ALL_ERRORS_FATAL TRUE)

project(Pets)

add_executable(pets graph/graph.cc 
                    pets/contract_graph.cc
                    pets/generate_pattern.cc 
                    pets/select_pattern.cc
                    pets/update_pattern.cc
                    pets/pets.cc)

target_include_directories(pets PRIVATE ${PROJECT_SOURCE_DIR})