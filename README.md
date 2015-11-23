This program is part of a research project about factoring integers via lattice algorithms.

###Compiling
You can compile this code via 

```
   cmake --build ./build --target FactoringIntegersViaLatticeAlgorithms -- -j 4
```

This program requires NTL >=9.5.0 and GMP to compile NTL thread safe. 
The link to NTL and GMP is hardcoded in the `CMakeLists.txt` file, so you
have to change that before you run cmake.

You can also change the number of threads, that will be used, in `CMakeLists.txt` by changing the line

```
add_definitions(-D__NUM_THREADS__=2)
```

By default the program will use two threads.

###Licence
The code is licenced under the GPLv2. See LICENCE.md for more details.