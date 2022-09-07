
add_executable(Navier_Stokes_steady_2D Navier_Stokes_steady_2D.cpp)
target_include_directories(Navier_Stokes_steady_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Navier_Stokes_steady_2D ${PYTHON_LIBRARIES})

add_executable(Burger Burger.cpp)
target_include_directories(Burger PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Burger ${PYTHON_LIBRARIES})
