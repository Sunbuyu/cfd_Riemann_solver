cmake_minimum_required(VERSION 3.20)
project(cfd_Riemann_solver)

set(CMAKE_CXX_STANDARD 14)

add_executable(cfd_Riemann_solver
        input.txt
        main.cpp)
