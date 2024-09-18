cmake_minimum_required(VERSION 3.29)
project(Call_1 C)

set(CMAKE_C_STANDARD 11)

find_package(Python3 COMPONENTS Interpreter Development NumPy)
set(PYTHON_INCLUDE_DIR /Users/mofefagade/anaconda3/include/python3.11)
include_directories(${PYTHON_INCLUDE_DIR})

add_executable(Call_1 main.c)
target_link_libraries(Call_1 ${Python3_LIBRARIES})