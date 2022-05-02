cmake_minimum_required(VERSION 3.21)
project(optfilmupdate)

set(CMAKE_CXX_STANDARD 17)

add_executable(Opt_FILMupdate main.cpp film.h pwlf.h data.h
        filmadastorage.h filmadalru.h
        )
