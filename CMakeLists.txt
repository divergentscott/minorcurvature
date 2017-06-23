cmake_minimum_required(VERSION 2.8)

PROJECT(lesscurve)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(lesscurve MACOSX_BUNDLE lesscurve.cxx)
#
  if(VTK_LIBRARIES)
    target_link_libraries(lesscurve ${VTK_LIBRARIES})
  else()
    target_link_libraries(lesscurve vtkHybrid vtkWidgets)
  endif()
