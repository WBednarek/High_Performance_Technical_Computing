cmake_minimum_required(VERSION 3.6)
project(hpc_assignment)

#set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(SOURCE_FILES

        hpc-assignment/hpc-assignment/explicit-upwind.cpp
        #hpc-assignment/hpc-assignment/implicit-upwind.cpp
        #hpc-assignment/hpc-assignment/crank-nicolson.cpp
        )

add_executable(hpc_assignment ${SOURCE_FILES})


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")


target_link_libraries(hpc_assignment ${MPI_LIBRARIES} m)

if (MPI_COMPILE_FLAGS)
    set_target_properties(hpc_assignment PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()
if (MPI_LINK_FLAGS)
    set_target_properties(hpc_assignment PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()