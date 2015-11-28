cmake_minimum_required(VERSION 3.1)
project(FactoringIntegersViaLatticeAlgorithms C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -lpthread -m64")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lpthread -m64")

set(SOURCE_FILES startProgram.cpp
        Equation.h Equation.cpp
        Factoring.h Factoring.cpp
        NewEnum.h NewEnum.cpp
        FileOutput.h FileOutput.cpp
        Statistics.h Statistics.cpp
        Timer.h Timer.cpp
        FactoringSettings.h OtherTests.h)

add_executable(FactoringIntegersViaLatticeAlgorithms ${SOURCE_FILES})

# Link NTL
set(NTL_LIB "C:/msys64/mingw64/lib/libntl.a")
target_link_libraries(FactoringIntegersViaLatticeAlgorithms ${NTL_LIB})

# Linker
set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc")

FIND_PACKAGE ( Threads REQUIRED )
TARGET_LINK_LIBRARIES ( FactoringIntegersViaLatticeAlgorithms ${CMAKE_THREAD_LIBS_INIT} )