cmake_minimum_required(VERSION 3.10)
project(VehicleMPC)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 编译选项
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")
else()
    set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra")
endif()

# 查找Eigen库
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "Found Eigen3: ${Eigen3_DIR}")
else()
    message(FATAL_ERROR "Eigen3 not found. Please install Eigen3 first.")
endif()

# 配置qpOASES库
# 假设qpOASES已安装在系统中，或者在thirdparty目录下
option(USE_SYSTEM_QPOASES "Use system-installed qpOASES" ON)

if(USE_SYSTEM_QPOASES)
    find_library(QPOASES_LIB qpOASES)
    find_path(QPOASES_INCLUDE_DIR qpOASES.hpp)
    
    if(QPOASES_LIB AND QPOASES_INCLUDE_DIR)
        message(STATUS "Found system qpOASES: ${QPOASES_LIB}")
    else()
        message(WARNING "System qpOASES not found, falling back to local version")
        set(USE_SYSTEM_QPOASES OFF)
    endif()
endif()

if(NOT USE_SYSTEM_QPOASES)
    # 假设qpOASES源代码位于thirdparty/qpOASES
    set(QPOASES_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/qpOASES")
    if(NOT EXISTS ${QPOASES_SRC_DIR})
        message(FATAL_ERROR "qpOASES source not found in ${QPOASES_SRC_DIR}. Please check the path.")
    endif()
    
    # qpOASES源文件列表
    file(GLOB QPOASES_SOURCES 
        "${QPOASES_SRC_DIR}/src/*.cpp"
        "${QPOASES_SRC_DIR}/src/BLASReplacement.cpp"
        "${QPOASES_SRC_DIR}/src/LAPACKReplacement.cpp"
    )
    
    # 编译qpOASES为静态库
    add_library(qpOASES STATIC ${QPOASES_SOURCES})
    target_include_directories(qpOASES PUBLIC "${QPOASES_SRC_DIR}/include")
    target_compile_definitions(qpOASES PUBLIC __NO_COPYRIGHT__)
    
    set(QPOASES_LIB qpOASES)
    set(QPOASES_INCLUDE_DIR "${QPOASES_SRC_DIR}/include")
endif()

# 主程序
add_executable(vehicle_mpc 
    vehicle_mpc.cpp
)

# 链接库
target_link_libraries(vehicle_mpc 
    PRIVATE 
    Eigen3::Eigen 
    ${QPOASES_LIB}
)

# 包含头文件目录
target_include_directories(vehicle_mpc
    PRIVATE 
    ${QPOASES_INCLUDE_DIR}
)

# 安装目标
install(TARGETS vehicle_mpc
    RUNTIME DESTINATION bin
)
    