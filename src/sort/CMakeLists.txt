include_directories(include)

file(GLOB SortImpl *.cpp)

# build this library as a static library
add_library(Sort STATIC ${SortImpl})