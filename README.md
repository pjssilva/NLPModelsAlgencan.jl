Algencan.jl
===========

**Algencan.jl** is a [Julia](http://julialang.org/) interface to the
[Algencan](https://www.ime.usp.br/~egbirgin/tango/codes.php)
nonlinear solver.

**Status**

At this point this is pre-alpha software. The main structure will change and
only unconstrained optimization is possible.

**Installation**:

I have not yet created a Package out of it. Up to now it is only a module,  so
to use it you just clone the repository and add the `src` directory to Julia's
`LOAD_PATH`, see the example and the `example` directory.

Algencan is a high performance and large scale Augmented Lagrangian solver
written by Ernesto Birgin and Mario Martínez. It has many special features like
being able to use HSL library functions to speed up linear algebra with sparse
matrices and some smart acceleration strategies. To use it, at least for now,
you need to compile Algencan yourself. Go to Tango's project
[website](https://www.ime.usp.br/~egbirgin/tango/codes.php) grab Algencan and
compile it. However, you need to make the small changes below to be able to get
a dynamic library from it as Algencan creates a static library by default.

1. Add the option `-f PIC` to  `CFLAGS` in the top of the main Makefile.

2. After compiling, you'll get a libalgencan.a file in the `lib` subdir. In
order to create a dynamic library try the following:

```
gcc -shared -o libalgencan.so -Wl,--whole-archive libalgencan.a \\
    -Wl,--no-whole-archive -lgfortran
```

It should create a file named `libalgencan.so` in the `lib` dir.

3. Create a environmental library named `ALGENCAN_LIB_DIR` pointing to the
`lib` dir path.
"""
