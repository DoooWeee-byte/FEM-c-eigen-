
add_executable(Stokes_Unsteady_2D Stokes_Unsteady_2D.cpp)
target_include_directories(Stokes_Unsteady_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Stokes_Unsteady_2D ${PYTHON_LIBRARIES})
