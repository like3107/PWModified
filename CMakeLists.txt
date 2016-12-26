cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(PowerWatershed)

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB src ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/include/*.h)
source_group(SOURCE FILES ${src})

file(GLOB pinkSrc ${PROJECT_SOURCE_DIR}/PINK/*.h ${PROJECT_SOURCE_DIR}/PINK/*.cpp)
source_group(PINK FILES ${pinkSrc})

file(GLOB csparse_src ${PROJECT_SOURCE_DIR}/CSparse/Source/*.cpp ${PROJECT_SOURCE_DIR}/CSparse/include/*.h)
list(REMOVE_ITEM csparse_src ${PROJECT_SOURCE_DIR}/CSparse/Source/cs_chol.cpp ${PROJECT_SOURCE_DIR}/CSparse/Source/cs_cholsol.cpp)
source_group(CSPARSE FILES ${csparse_src})
include_directories(${PROJECT_SOURCE_DIR}/CSparse)
include_directories(${PROJECT_SOURCE_DIR}/CSparse/include)

file(GLOB argv_src ${PROJECT_SOURCE_DIR}/argv/*.h ${PROJECT_SOURCE_DIR}/argv/*.cpp)
source_group(ARGV FILES ${argv_src})

#opencv
find_package(OpenCV REQUIRED)

#itk
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#glog
include("${PROJECT_SOURCE_DIR}/cmake/External/glog.cmake")
include_directories(${GLOG_INCLUDE_DIRS})

#gflags
include("${PROJECT_SOURCE_DIR}/cmake/External/gflags.cmake")
include_directories(${GFLAGS_INCLUDE_DIRS})

add_executable(PowerWatershed ${src} ${pinkSrc} ${csparse_src} ${argv_src})
target_link_libraries(PowerWatershed ${OpenCV_LIBS} ${ITK_LIBRARIES} ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})


