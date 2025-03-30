cmake_minimum_required(VERSION 3.10)

# set the compiler flags
add_compile_options("-Wall")

# set up the project itself
project(RadixSort LANGUAGES CXX)

# adding the src subdirectory, as we have code to build there
add_subdirectory(src)
