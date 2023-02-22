

[![DOI](https://zenodo.org/badge/30722436.svg)](https://zenodo.org/badge/latestdoi/30722436)



This program is part of a research project about factoring integers via lattice algorithms.

### Compiling
You can compile this code via

```
cmake --build ./build --target FactoringIntegersViaLatticeAlgorithms -- -j 4
```

This program uses NTL. The link to the library is hardcoded in the `CMakeLists.txt` file, so you
have to change that before you run cmake. Maybe you also have to add
`include_directories(path_to_ntl_include)` to `CMakeLists.txt`.

### Results
The directory `results` contains a database dump from a MySQL 5.6 database containing all
results of the 1680 tests done during the project. Have a look at QUERIES.md for the database
structure and some example queries.
The files generated during the tests are not included due to their size. They sum up to about 2 GB
and even compressed they need more than 600 MB, so they can not be hosted on GitHub.

### Licence
The code is licenced under the GPLv2. See LICENCE.md for more details.

The results in the results directory are licenced under the Creative-Commons Licence CC-BY-NC-SA 4.0. Have a look at http://creativecommons.org/licenses/by-nc-sa/4.0/ for details.
