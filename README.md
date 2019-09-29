# divide-and-conquer-eigenvalues

This is realization of divide-and-conquer eigenvalues algorithm for symmetric tridiagonal matrix, designed by Cuppen in 1982.

### Tech

Algorithm is coded in pure C11.
It is single-threaded, but can be simply parallelized, giving recursive tasks for T1 and T2 matricies to thread-pool.

### Installation

`divide-and-conquer-eigenvalues` requires only `C89`-compatible compiler and `make` utility.

```sh
$ cd divide-and-conquer-eigenvalues
$ make
$ ./bin/test data/in.mat
```
