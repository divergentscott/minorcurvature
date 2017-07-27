cmake_minimum_required(VERSION 2.8)
set (CMAKE_CXX_STANDARD 11)
include(ExternalProject)

PROJECT(lesscurve)

find_package(VTK REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include(${VTK_USE_FILE})

set(BLOSSOM_DIR blossom5-v2.05.src)
file(GLOB_RECURSE BLOSSOM_SRC ${BLOSSOM_DIR}/*.h  ${BLOSSOM_DIR}/*.cpp)


add_executable(lesscurve MACOSX_BUNDLE lesscurve.cxx)

add_executable(graph_writing MACOSX_BUNDLE graph_writing.cxx ${BLOSSOM_SRC})


  if(VTK_LIBRARIES)
    target_link_libraries(lesscurve ${VTK_LIBRARIES})
    target_link_libraries(graph_writing ${VTK_LIBRARIES})
  else()
    target_link_libraries(lesscurve vtkHybrid vtkWidgets)
    target_link_libraries(graph_writing vtkHybrid vtkWidgets)
  endif()

target_link_libraries (lesscurve Eigen3::Eigen)
target_link_libraries (graph_writing Eigen3::Eigen)
