cmake_minimum_required(VERSION 3.8)
project(test_ikfom)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
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

if(BUILD_TESTING)
  add_executable(test_ikfom src/test_ikfom.cpp)
  target_link_libraries(test_ikfom Eigen3::Eigen)

  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
