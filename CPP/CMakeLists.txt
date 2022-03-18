cmake_minimum_required(VERSION 2.8)
project(main)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(sources_dir "${CMAKE_CURRENT_SOURCE_DIR}")
# set(sources "${sources_dir}/src/main.cpp ${sources_dir}/src/ReactionDiffusion.cpp ${sources_dir}/src/ForwardEuler.cpp")
set(sources "${sources_dir}/src/main.cpp")

add_executable(main ${sources})

find_package(BLAS REQUIRED)
find_package(Boost REQUIRED)

target_link_libraries(main ${BLAS_LIBRARIES})
target_link_libraries(main ${Boost_LIBRARIES})