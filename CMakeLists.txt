cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)
PROJECT(lesscurve)
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) 
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) 
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin) 

include(conan)
conan_cmake_run(CONANFILE conanfile.py
		BASIC_SETUP CMAKE_TARGETS 
		BUILD missing
		)
		
set(BLOSSOM_DIR blossom5-v2.05.src)
file(GLOB_RECURSE BLOSSOM_SRC ${BLOSSOM_DIR}/*.h  ${BLOSSOM_DIR}/*.cpp)

# include_directories(${CMAKE_SOURCE_DIR})

add_executable(lesscurve lesscurve.cxx ${BLOSSOM_SRC})
target_link_libraries(lesscurve CONAN_PKG::VTK-maxi CONAN_PKG::eigen CONAN_PKG::boost)
#add_executable(lesscurve lesscurve.cxx ${BLOSSOM_SRC})

# add_executable(graph_writing graph_writing.cxx ${BLOSSOM_SRC})
# target_link_libraries(graph_writing CONAN_PKG::VTK-maxi CONAN_PKG::eigen CONAN_PKG::boost)
# target_link_libraries(lesscurve ${CMAKE_DL_LIBS} CONAN_PKG::VTK CONAN_PKG::boost CONAN_PKG::eigen)