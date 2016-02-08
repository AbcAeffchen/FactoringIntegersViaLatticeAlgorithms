cmake_minimum_required(VERSION 3.1)
project(FactoringIntegersViaLatticeAlgorithms C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -lpthread -m64")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lpthread -m64")

set(SOURCE_PATH src)
set(SOURCE_FILES ${SOURCE_PATH}/startProgram.cpp
        ${SOURCE_PATH}/Equation.h ${SOURCE_PATH}/Equation.cpp
        ${SOURCE_PATH}/Factoring.h ${SOURCE_PATH}/Factoring.cpp
        ${SOURCE_PATH}/NewEnum.h ${SOURCE_PATH}/NewEnum.cpp
        ${SOURCE_PATH}/FileOutput.h ${SOURCE_PATH}/FileOutput.cpp
        ${SOURCE_PATH}/Statistics.h ${SOURCE_PATH}/Statistics.cpp
        ${SOURCE_PATH}/Timer.h ${SOURCE_PATH}/Timer.cpp
        ${SOURCE_PATH}/FactoringSettings.h ${SOURCE_PATH}/OtherTests.h)

add_executable(FactoringIntegersViaLatticeAlgorithms ${SOURCE_FILES})

#use centered continued fractions
#add_definitions(-DFS_CCF)

# Link NTL
set(NTL_LIB "C:/msys64/mingw64/lib/libntl.a")
target_link_libraries(FactoringIntegersViaLatticeAlgorithms ${NTL_LIB})

# Linker
set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc")

FIND_PACKAGE ( Threads REQUIRED )
TARGET_LINK_LIBRARIES ( FactoringIntegersViaLatticeAlgorithms ${CMAKE_THREAD_LIBS_INIT} )