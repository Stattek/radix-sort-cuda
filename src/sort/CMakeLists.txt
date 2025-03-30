include_directories(include)

file(GLOB SortImpl *.cpp)

add_library(Sort STATIC ${SortImpl})