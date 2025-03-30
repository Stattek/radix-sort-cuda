add_subdirectory(sort)

file(GLOB SrcImpl *.cpp)

add_executable(RadixSort ${SrcImpl})

target_link_libraries(Sort)
