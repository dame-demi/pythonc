cmake_minimum_required(VERSION 3.22)
project(untitled2 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

include_directories(include)

add_executable(untitled2 main.c)

# Set the path to the SUNDIALS installation directory
set(SUNDIALS_DIR "/home/tony/Downloads/instdir/lib/cmake/sundials")

# Using Sundials as a Third Party Library
find_package(SUNDIALS 7.1.0 REQUIRED)

# Find the math library
find_library(MATH_LIBRARY m)

# Link the SUNDIALS libraries
target_link_libraries(untitled2 PUBLIC SUNDIALS::nvecserial SUNDIALS::core SUNDIALS::sunmatrixdense SUNDIALS::sunmatrixband SUNDIALS::sunlinsoldense SUNDIALS::sunlinsolband SUNDIALS::cvode ${MATH_LIBRARY})