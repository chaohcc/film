cmake_minimum_required(VERSION 3.21)
project(FILM)

set(CMAKE_CXX_STANDARD 17)

add_executable(FILM main.cpp film.h pwlf.h data.h
        filmadastorage.h filmadalru.h
        )

if(MSVC)
    set(CMAKE_CXX_FLAGS "/O2 /arch:AVX2 /W1 /EHsc")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "-O3 -xHost")
else()
    # clang and gcc
    set(CMAKE_CXX_FLAGS "-O0 -march=native -Wall -Wextra")    #chao changes O3 to O0
endif()
