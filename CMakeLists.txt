cmake_minimum_required(VERSION 3.1)
project(FactoringIntegersViaLatticeAlgorithms C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3  -lpthread")

set(SOURCE_FILES startProgram.cpp
                 Equation.h Equation.cpp
                 Factoring.h Factoring.cpp
                 FileOutput.h FileOutput.cpp
                 NewEnum.h NewEnum.cpp
                 Statistics.h Statistics.cpp
                 Timer.h Timer.cpp)

add_executable(FactoringIntegersViaLatticeAlgorithms ${SOURCE_FILES})

# Link NTL
set(NTL_LIB "C:/Program Files/mingw-w64/x86_64-4.9.1-posix-seh-rt_v3-rev0/mingw64/x86_64-w64-mingw32/lib/libntl.a")
target_link_libraries(FactoringIntegersViaLatticeAlgorithms ${NTL_LIB})

# Linker
set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc")

FIND_PACKAGE ( Threads REQUIRED )
TARGET_LINK_LIBRARIES ( FactoringIntegersViaLatticeAlgorithms ${CMAKE_THREAD_LIBS_INIT} )