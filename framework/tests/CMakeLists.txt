cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES test_main.cpp kernel_tests.cpp physics_tests.cpp particle_tests.cpp)
    
add_executable(learnSPH_tests ${SOURCE_FILES})
add_dependencies(learnSPH_tests learnSPH)
target_link_libraries(learnSPH_tests learnSPH)
