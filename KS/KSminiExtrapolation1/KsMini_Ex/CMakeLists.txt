
add_executable(KsMini_Ex KsMini_Ex.cpp)
target_include_directories(KsMini_Ex PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(KsMini_Ex ${PYTHON_LIBRARIES})


add_executable(KsMini_Ex_Steady KsMini_Ex_Steady.cpp)
target_include_directories(KsMini_Ex_Steady PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(KsMini_Ex_Steady ${PYTHON_LIBRARIES})
