add_executable(AllenChan_2D AllenChan_2D.cpp)
target_include_directories(AllenChan_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(AllenChan_2D ${PYTHON_LIBRARIES})


