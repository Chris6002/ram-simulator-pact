cmake_minimum_required(VERSION 2.8)
project(ram_simulator)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(OpenCV REQUIRED)
add_executable(ram_simulator ram_simulator.cpp)
target_link_libraries(ram_simulator ${OpenCV_LIBS})
