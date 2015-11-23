This program is part of a research project about factoring integers via lattice algorithms.

###Compiling
You can compile this code via 

```
   cmake --build ./build --target FactoringIntegersViaLatticeAlgorithms -- -j 4
```

This program uses NTL. The link to the library is hardcoded in the `CMakeLists.txt` file, so you
have to change that before you run cmake. Maybe you also have to add 
`include_directories(path_to_ntl_include)` to `CMakeLists.txt`.

###Licence
The code is licenced under the GPLv2. See LICENCE.md for more details.