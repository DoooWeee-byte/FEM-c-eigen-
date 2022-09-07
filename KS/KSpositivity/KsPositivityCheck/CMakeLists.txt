

#add_executable(KsPosCheck KsPosCheck.cpp)
#target_include_directories(KsPosCheck PRIVATE ${PYTHON_INCLUDE_DIRS})
#target_link_libraries(KsPosCheck ${PYTHON_LIBRARIES})

add_executable(Kserrorcheck_11 Kserrorcheck_11.cpp)
target_include_directories(Kserrorcheck_11 PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Kserrorcheck_11 ${PYTHON_LIBRARIES})



add_executable(Kserrorcheck1_plt Kserrorcheck1_plt.cpp)
target_include_directories(Kserrorcheck1_plt PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(Kserrorcheck1_plt ${PYTHON_LIBRARIES})

add_executable(plotOrder plotOrder.cpp)
target_include_directories(plotOrder PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(plotOrder ${PYTHON_LIBRARIES})
