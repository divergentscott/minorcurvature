cmake_minimum_required(VERSION 2.8)
set (CMAKE_CXX_STANDARD 11)

PROJECT(lesscurve)

find_package(VTK REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include(${VTK_USE_FILE})


add_executable(lesscurve MACOSX_BUNDLE lesscurve.cxx)
#
  if(VTK_LIBRARIES)
    target_link_libraries(lesscurve ${VTK_LIBRARIES})
  else()
    target_link_libraries(lesscurve vtkHybrid vtkWidgets)
  endif()
target_link_libraries (lesscurve Eigen3::Eigen)
