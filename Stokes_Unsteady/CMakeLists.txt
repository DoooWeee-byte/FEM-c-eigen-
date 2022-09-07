cmake_minimum_required(VERSION 3.0)

project(PDWFE)
include_directories(${PROJECT_SOURCE_DIR}/include)
message(STATUS "${PROJECT_SOURCE_DIR}")

find_package(PythonLibs 3)

message(STATUS "${PYHON_INCLUDE_DIRS}")

add_subdirectory(Stokes_Unsteady_2D)
