# cmake_minimum_required(VERSION 3.17)
# project(helibTest)  # Replace "YourProjectName" with your actual project name

# set(CMAKE_CXX_STANDARD 17)  # Or use the appropriate C++ version
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# # Add all your source files here
# add_executable(main main.cpp comparator.cpp tools.cpp common.cpp common.h host_data.h host_data_io.cpp host_data_io.h host_kernel.cpp host_kernel.h)  # Add all your .cpp files here

# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")

# # Set the path to the TFHE library
# set(GMP_LIBRARY "/home/faiz007/Desktop/FYP/HElib/helib_install/helib_pack/")
# set(NTL_LIBRARY "/home/faiz007/Desktop/FYP/HElib/helib_install/helib_pack/")
# set(GMP_LIBRARY_DIR "/home/faiz007/Desktop/FYP/HElib/helib_install/helib_pack/")
# set(HELIB_LIBRARY_DIR "/home/faiz007/Desktop/FYP/HElib/helib_install/helib_pack/")
# set(MINIMAP_LIBRARY "/home/faiz007/Desktop/FYP/HElib/helib_install/helib_pack/")

# # Add the TFHE include directory to the header search paths
# include_directories(${HELIB_LIBRARY_DIR}/include)
# # Add the directory containing minimap.h to the header search paths
# include_directories(${MINIMAP_LIBRARY}/lib/minimap2)

# # Link your executable with the TFHE library
# target_link_libraries(main
#         "${HELIB_LIBRARY_DIR}/lib/libhelib.a"
#         "${GMP_LIBRARY_DIR}/lib/libgmp.so.10"
#         "${GMP_LIBRARY_DIR}/lib/libntl.so.44"
#         "${GMP_LIBRARY}/lib/libgmp.so"
#         "${NTL_LIBRARY}/lib/libntl.so"
#         "${MINIMAP_LIBRARY}lib/minimap2/minimap.h"
# )

cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(comparison_circuit
       	LANGUAGES CXX)

find_package(helib REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})

add_executable(main main.cpp comparator.cpp tools.cpp)
add_executable(bgv_fatboot bgv_fatboot.cpp)

# add_executable(comparison_circuit comparison_circuit.cpp comparator.cpp tools.cpp)
# add_executable(sorting_circuit sorting_circuit.cpp comparator.cpp tools.cpp)
# add_executable(min_max_circuit min_max_circuit.cpp comparator.cpp tools.cpp)

target_link_libraries(main helib)
target_link_libraries(bgv_fatboot helib)
# target_link_libraries(bgv_fatboot benchmark::benchmark)
# target_link_libraries(comparison_circuit helib)
# target_link_libraries(sorting_circuit helib)
# target_link_libraries(min_max_circuit helib)

