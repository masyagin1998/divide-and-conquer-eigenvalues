# divide-and-conquer-eigenvalues

This is realization of divide-and-conquer eigenvalues algorithm for symmetric tridiagonal matrix, designed by Cuppen in 1982.

### Tech

Algorithm is coded in pure `C99`. You can compile it with `C89` compiler easily, if You change `long double` to `doule`, `sqrtl` and `fabsl` to `sqrt` and `abs`, or declare Your own `sqrtl` and `fabsl` (`C89` doesn't have them) and change C version in Makefile.

### Installation

`divide-and-conquer-eigenvalues` requires only `C99`-compatible compiler and `make` utility.

```sh
$ cd divide-and-conquer-eigenvalues
$ make
$ ./bin/test data/in.mat
```

### Short description

`divide-and-conquer-eigenvalues` finds for given symmetric tridiagonal matrix `T` matrices of eigenvectors `Q` and eigenvalues `L`, such `T = Q * L * Q_t`, partitioning large `T` matrix into two half-sized matrices `T1` and `T2` (they are also partitioned, so we use `divide-and-conquer` strategy).

### Optimization

You can add thread-pool to this algorithm:
- recursive `divide-and-conquer` steps can be processed in parallel. You can divide recursive calls to current thread and additional thread;
- `secular equation` can be solved in parallel too, because its routine doesn't modify any input matrices;

### Usefull articles

In english:
- [Arbenz P. Lecture Notes on Solving Large Scale Eigenvalue Problems](https://yadi.sk/i/7Ry-GgZ_vYoDvw)
- [Arbenz P. - his website with lectures, presentations, etc](http://people.inf.ethz.ch/arbenz/ewp/index.html)
- [Demmel J. W. Applied Numerical Linear Algebra](https://yadi.sk/i/CAk-jcKLxsj85Q)
- [Rutter J. A Serial Implementation of Cuppen's Divide and Conquer Algorithm for the Symmetric Eigenvalue Problem](https://yadi.sk/i/_UFPGDnKvpdfqA)

In russian:
- [Golubkov A. Yu. Lecture notes on computational methods of linear algebra.](https://yadi.sk/i/SQ0BW2aAu0m-LQ)
- [AlgoWiki Divide-and-Conquer](https://algowiki-project.org/ru/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%C2%AB%D1%80%D0%B0%D0%B7%D0%B4%D0%B5%D0%BB%D1%8F%D0%B9_%D0%B8_%D0%B2%D0%BB%D0%B0%D1%81%D1%82%D0%B2%D1%83%D0%B9%C2%BB_%D0%B2%D1%8B%D1%87%D0%B8%D1%81%D0%BB%D0%B5%D0%BD%D0%B8%D1%8F_%D1%81%D0%BE%D0%B1%D1%81%D1%82%D0%B2%D0%B5%D0%BD%D0%BD%D1%8B%D1%85_%D0%B7%D0%BD%D0%B0%D1%87%D0%B5%D0%BD%D0%B8%D0%B9_%D0%B8_%D0%B2%D0%B5%D0%BA%D1%82%D0%BE%D1%80%D0%BE%D0%B2_%D1%81%D0%B8%D0%BC%D0%BC%D0%B5%D1%82%D1%80%D0%B8%D1%87%D0%BD%D0%BE%D0%B9_%D1%82%D1%80%D0%B5%D1%85%D0%B4%D0%B8%D0%B0%D0%B3%D0%BE%D0%BD%D0%B0%D0%BB%D1%8C%D0%BD%D0%BE%D0%B9_%D0%BC%D0%B0%D1%82%D1%80%D0%B8%D1%86%D1%8B)
