

add_executable(FE_elasticity FE_elasticity.cpp)
target_include_directories(FE_elasticity PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(FE_elasticity ${PYTHON_LIBRARIES})
