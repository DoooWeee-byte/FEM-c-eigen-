add_executable(FE_heat_2D FE_heat_2D.cpp)
target_include_directories(FE_heat_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_heat_2D ${PYTHON_LIBRARIES})




add_executable(FE_wave_2D FE_wave_2D.cpp)
target_include_directories(FE_wave_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_wave_2D ${PYTHON_LIBRARIES})
