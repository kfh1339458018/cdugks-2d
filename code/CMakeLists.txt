# cmake minimum required
cmake_minimum_required(VERSION 3.10)

# define project name and its language
project(cdugks-2d LANGUAGES CXX) 

# set MPICH for C and C++
set(MPI_C_COMPILER /home/npu/Software_Install/mpich-3.2.1/bin/mpicc)
set(MPI_CXX_COMPILER /home/npu/Software_Install/mpich-3.2.1/bin/mpicxx)

# find MPI package
find_package(MPI REQUIRED)

# set C++ standard
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# add executable file
add_executable(cdugks_2d
    main.cpp 
    cdugks_2d.cpp 
    config.cpp 
    mesh_read.cpp 
    output.cpp 
    mpi_process.cpp 
    test_cylinder_half.cpp 
    test_cylinder.cpp 
    test_cavity.cpp 
    test_riemann_2d.cpp 
    test_riemann_1d.cpp 
)


include_directories(${MPI_INCLUDE_PATH})


target_link_libraries(cdugks_2d ${MPI_CXX_LIBRARIES})

