# compile this subdirectory
add_subdirectory(sort)

# glob all the files to compile
file(GLOB SrcImpl *.cpp)

# build the executable
add_executable(RadixSort ${SrcImpl})

# link this executable to the Sort library in the subdirectory
target_link_libraries(Sort)
