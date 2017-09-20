cmake_minimum_required(VERSION 3.1)
project(FactoringIntegersViaLatticeAlgorithms C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -lpthread -m64")

set(SOURCE_PATH src)
set(SOURCE_FILES ${SOURCE_PATH}/startProgram.cpp
        ${SOURCE_PATH}/FileOutput.h ${SOURCE_PATH}/FileOutput.cpp
        ${SOURCE_PATH}/FactoringSettings.h ${SOURCE_PATH}/OtherTests.h)

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