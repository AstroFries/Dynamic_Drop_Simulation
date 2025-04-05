cmake_minimum_required(VERSION 3.14)

project(
  "DynamicDrop"
  VERSION 1.0
  LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

message(STATUS "Started CMake for ${PROJECT_NAME} v${PROJECT_VERSION}...\n")
#set(BUILD_SHARED_LIBS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

set(FRAMEWORK2D_DIR "${PROJECT_SOURCE_DIR}")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(THIRD_PARTY_DIR "${PROJECT_SOURCE_DIR}/third_party")
set(BINARY_DIR "${PROJECT_SOURCE_DIR}/bin")
set(LIBRARY_DIR "${PROJECT_SOURCE_DIR}/libs")

add_subdirectory(src)
add_subdirectory(third_party)

# 获取src目录下的所有源文件
file(GLOB_RECURSE SOURCES "src/*.cpp")
list(FILTER SOURCES EXCLUDE REGEX "src/Reader/.*")

file(GLOB_RECURSE READER_SOURCES "src/Reader/*.cpp")

# 添加可执行文件
add_executable(${PROJECT_NAME} ${SOURCES})
add_executable(Reader ${READER_SOURCES})

#Eigen
set(EIGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/Eigen")
target_include_directories(DynamicDrop PRIVATE ${EIGEN_INCLUDE_DIR})

set(Python3_ROOT_DIR "C:/Users/12989/AppData/Local/Programs/Python/Python310")
set(Python3_INCLUDE_DIRS "C:/Users/12989/AppData/Local/Programs/Python/Python310/include")
set(Python3_LIBRARIES "C:/Users/12989/AppData/Local/Programs/Python/Python310/libs/python310.lib")
#set(Python_ROOT_DIR "C:/Users/12989/AppData/Local/Programs/Python/Python310")

target_include_directories(DynamicDrop PRIVATE
    ${Python3_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/third_party
)

target_link_libraries(DynamicDrop PRIVATE
    ${Python3_LIBRARIES}  # 链接python310.lib等库
)

set(NUM_PY_INCLUDE_DIR "C:/Users/12989/AppData/Local/Programs/Python/Python310/Lib/site-packages/numpy/_core/include")
target_include_directories(DynamicDrop PRIVATE ${NUM_PY_INCLUDE_DIR})

message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")
set(SFML_DIR "external/SFML-3.0.0/lib/cmake/SFML")
find_package(SFML 3 REQUIRED COMPONENTS Graphics Audio Network)
target_link_libraries(DynamicDrop PRIVATE SFML::Graphics SFML::Audio SFML::Network SFML::Window)
target_link_libraries(Reader PRIVATE SFML::Graphics SFML::Audio SFML::Network SFML::Window)

#data
set(DATA_DIR "${PROJECT_SOURCE_DIR}/data")
target_compile_definitions(${PROJECT_NAME} PRIVATE DATA_DIR="${DATA_DIR}")
target_compile_definitions(Reader PRIVATE DATA_DIR="${DATA_DIR}")