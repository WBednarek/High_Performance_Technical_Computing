cmake_minimum_required(VERSION 3.5)
project(High_Performance_Technical_Computing)

#set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

#[[
Usefull flags: -Wall -pedantic -Wextra -Weffc++ the last checking if code is written accorting to books:
"Effective C++" and "More Effective C++" written by Scott Meyers
]]

set(SOURCE_FILES
        Display.cpp
        Display.h
        ExplicitUpwindScheme.cpp
        ExplicitUpwindScheme.h
        GeneralScheme.cpp
        GeneralScheme.h
        ImplicitUpwindScheme.cpp
        ImplicitUpwindScheme.h
        Lax_Wendroff.cpp
        Lax_Wendroff.h
        main.cpp
        main.h
        MathFunctions.cpp
        MathFunctions.h
        Matrix.cpp
        Matrix.h
        Richtmyer_multi_step.cpp
        Richtmyer_multi_step.h
        ExplicitUpwindParallel.cpp
        ExplicitUpwindParallel.h
        ImplicitParallel.h
        ImplicitParallel.cpp ExplicitUpwindParallel.cpp ExplicitUpwindParallel.h CrankParallel.cpp CrankParallel.h)

add_executable(High_Performance_Technical_Computing ${SOURCE_FILES})


target_link_libraries(High_Performance_Technical_Computing ${MPI_LIBRARIES} m)

if (MPI_COMPILE_FLAGS)
    set_target_properties(High_Performance_Technical_Computing PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()
if (MPI_LINK_FLAGS)
    set_target_properties(High_Performance_Technical_Computing PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()