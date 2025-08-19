cmake_minimum_required(VERSION 3.15)
project(PortfolioRiskAnalyzer CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Eigen (header-only). Option A: system package. Option B: vendor as subdir.
# If installed system-wide, something like:
# find_package(Eigen3 3.4 REQUIRED NO_MODULE)
# target_link_libraries(risk PRIVATE Eigen3::Eigen)
# Otherwise, vendor Eigen in third_party/eigen and add_subdirectory.

add_executable(risk main.cpp)
target_include_directories(risk PRIVATE ${CMAKE_SOURCE_DIR})
