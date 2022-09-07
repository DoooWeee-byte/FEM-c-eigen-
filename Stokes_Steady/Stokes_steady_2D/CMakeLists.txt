
add_executable(Stokes_steady_2D Stokes_steady_2D.cpp)
target_include_directories(Stokes_steady_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Stokes_steady_2D ${PYTHON_LIBRARIES})


add_executable(Stokes_steady_2D_model1 Stokes_steady_2D_model1.cpp)
target_include_directories(Stokes_steady_2D_model1 PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Stokes_steady_2D_model1 ${PYTHON_LIBRARIES})
