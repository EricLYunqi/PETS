cmake_minimum_required(VERSION 3.9)
project(TATTOO)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -pthread")

find_package(
        Boost
        REQUIRED
        COMPONENTS graph thread
)

set(UTILS_SOURCE_FILES
        utils/json.hpp
        utils/NetSimile.h
        utils/NetSimile.cpp)

set(TATTOO_SOURCE_FILES
        tattoo.cpp
        FastPGGraph.h
        FastPGGraph.cpp
        base.h
        base.cpp
        ${UTILS_SOURCE_FILES})
add_executable(tattoo ${TATTOO_SOURCE_FILES})
target_link_libraries(tattoo ${Boost_LIBRARIES})

set(DAVINCI_SOURCE_FILES
        Davinci.cpp
        FastPGGraph.h
        FastPGGraph.cpp
        base.h
        base.cpp
        ${UTILS_SOURCE_FILES})
add_executable(Davinci ${DAVINCI_SOURCE_FILES})
target_link_libraries(Davinci ${Boost_LIBRARIES})
