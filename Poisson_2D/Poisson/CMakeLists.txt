add_executable(FE_Dirichlet_2D FE_Dirichlet_2D.cpp)
target_include_directories(FE_Dirichlet_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_Dirichlet_2D ${PYTHON_LIBRARIES})


add_executable(FE_Di_Neumann_2D FE_Di_Neumann_2D.cpp)
target_include_directories(FE_Di_Neumann_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_Di_Neumann_2D ${PYTHON_LIBRARIES})


add_executable(FE_Di_Robin_2D FE_Di_Robin_2D.cpp)
target_include_directories(FE_Di_Robin_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_Di_Robin_2D ${PYTHON_LIBRARIES})
