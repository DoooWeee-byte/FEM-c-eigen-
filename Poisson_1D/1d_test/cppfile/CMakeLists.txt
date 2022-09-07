add_executable(Poisson1D Poisson1D.cpp)
target_include_directories(Poisson1D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Poisson1D ${PYTHON_LIBRARIES})
