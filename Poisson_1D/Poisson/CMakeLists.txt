add_executable(FE_Dirich_Neu_1D FE_Dirich_Neu_1D.cpp)
target_include_directories(FE_Dirich_Neu_1D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_Dirich_Neu_1D ${PYTHON_LIBRARIES})
