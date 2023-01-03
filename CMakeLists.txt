###############################################
# CMake options
###############################################

cmake_minimum_required (VERSION 3.5.1 FATAL_ERROR)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(GL_vs_VK)

###############################################
# GLvsVK
###############################################

add_subdirectory(project)
