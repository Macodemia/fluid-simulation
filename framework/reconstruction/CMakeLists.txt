cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES main.cpp)
   
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

add_executable(learnSPH_reconstruction ${SOURCE_FILES})
add_dependencies(learnSPH_reconstruction learnSPH)
target_link_libraries(learnSPH_reconstruction learnSPH)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(learnSPH_reconstruction Threads::Threads)
