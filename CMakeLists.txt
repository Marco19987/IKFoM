cmake_minimum_required(VERSION 3.8)
project(ikfom_toolkit)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()


# find dependencies
find_package(ament_cmake REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)

include_directories(include)
ament_export_include_directories("include")

ament_export_dependencies(Eigen3)



install(
  DIRECTORY include/
  DESTINATION include
)

add_executable(test_ikfom test_ikfom/test_ikfom.cpp)
target_link_libraries(test_ikfom Eigen3::Eigen)

add_executable(test_ikfom_dyn_measurement test_ikfom/test_ikfom_dyn_measurement.cpp)
target_link_libraries(test_ikfom_dyn_measurement Eigen3::Eigen)


ament_package()
