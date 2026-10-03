cmake_minimum_required(VERSION 3.15)
project(TutorialAlgebra VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(AriadneAlgebra REQUIRED)

add_executable(tutorial_algebra tutorial_algebra.cpp)
target_link_libraries(tutorial_algebra PRIVATE AriadneAlgebra::ariadne-algebra)
