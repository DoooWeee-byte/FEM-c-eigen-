add_executable(KS_2D_test1 KS_2D_test1.cpp)
target_include_directories(KS_2D_test1 PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(KS_2D_test1 ${PYTHON_LIBRARIES})


add_executable(FE_Dirichlet_2D FE_Dirichlet_2D.cpp)
target_include_directories(FE_Dirichlet_2D PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_Dirichlet_2D ${PYTHON_LIBRARIES})



add_executable(Eigen_test1 Eigen_test1.cpp)
target_include_directories(Eigen_test1 PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Eigen_test1 ${PYTHON_LIBRARIES})
