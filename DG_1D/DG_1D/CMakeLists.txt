add_executable(DG_1D DG_1D.cpp)
target_include_directories(DG_1D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(DG_1D ${PYTHON_LIBRARIES})
