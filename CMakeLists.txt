cmake_minimum_required(VERSION 3.29)
project(Python_C C)

set(CMAKE_C_STANDARD 11)

find_package(Python3 COMPONENTS Interpreter Development NumPy)
set(PYTHON_INCLUDE_DIR /Users/mofefagade/anaconda3/include/python3.11)
include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${Python3_INCLUDE_DIR})

add_executable(Python_C /src/main.c /src/mytest.py)
target_link_libraries(Python_C ${Python3_LIBRARIES})