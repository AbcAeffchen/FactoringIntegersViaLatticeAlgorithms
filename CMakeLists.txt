cmake_minimum_required(VERSION 3.1)
project(FactoringIntegersViaLatticeAlgorithms C CXX)

set(CMAKE_CONFIGURATION_TYPES "EXPERIMENT;DEBUG" CACHE STRING "" FORCE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -m64 -fpermissive -Wall -Wextra -Wold-style-cast -Wuseless-cast")

if(DEFINED CMAKE_BUILD_TYPE)
else()
    set(CMAKE_BUILD_TYPE "EXPERIMENT")
endif()
message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE MATCHES "EXPERIMENT")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
endif()

set(SOURCE_FILES src/startProgram.cpp)

add_executable(FactoringIntegersViaLatticeAlgorithms ${SOURCE_FILES})

#use centered continued fractions
#add_definitions(-DFS_CCF)

# Link NTL
#set(NTL_LIB "C:/msys64/mingw64/lib/libntl.a")
#target_link_libraries(FactoringIntegersViaLatticeAlgorithms ${NTL_LIB})
target_link_libraries(FactoringIntegersViaLatticeAlgorithms ntl m gmp)

# Linker
#set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc")

FIND_PACKAGE ( Threads REQUIRED )
TARGET_LINK_LIBRARIES ( FactoringIntegersViaLatticeAlgorithms ${CMAKE_THREAD_LIBS_INIT} )