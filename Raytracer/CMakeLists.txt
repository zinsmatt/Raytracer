project(Raytracer)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)

find_package(SDL REQUIRED)

set(HEADER raytracer.h)


add_executable(${PROJECT_NAME} ${SRC_LIST} ${HEADER})

target_link_libraries(${PROJECT_NAME} SDL ${SDL_LIBRARIES})

set(CMAKE_CXX_FLAGS "-std=c++11 -O0")
