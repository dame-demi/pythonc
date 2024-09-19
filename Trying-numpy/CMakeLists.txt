cmake_minimum_required(VERSION 3.29)
project(Trying_numpy C)

set(CMAKE_C_STANDARD 11)

find_package(Python3 COMPONENTS Interpreter Development NumPy)
set(PYTHON_INCLUDE_DIR /Users/mofefagade/anaconda3/include/python3.11)
include_directories(${PYTHON_INCLUDE_DIR})

add_executable(Trying_numpy main_2.c)
target_link_libraries(Trying_numpy ${Python3_LIBRARIES})