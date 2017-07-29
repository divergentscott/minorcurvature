cmake_minimum_required(VERSION 2.8)
# set(CMAKE_C_COMPILER /usr/bin/clang)
# set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_BUILD_TYPE Release)
# # find_package(OpenMP)
# # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
# # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKR_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
# # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#
# set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
#
#
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
