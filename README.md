# Algencan.jl

**Algencan.jl** is a [JuMP / MathProgBase](https://www.juliaopt.org/) interface
to the [Algencan](https://www.ime.usp.br/~egbirgin/tango/codes.php)
nonlinear solver.

Algencan is a high performance and large scale Augmented Lagrangian solver
written by Ernesto Birgin and Mario Martínez. It has many special features like
being able to use HSL library functions to speed up linear algebra with sparse
matrices and some smart acceleration strategies.

## Status

At this point this is alpha software. From verision v0.2.0 on the code will work
with Julia 0.7 or later only. Use Julia 1.0  or later, if possible. If you need
to run Algencan.jl with the old Julia 0.6, please install the version v0.1.x.

## Installation

It is now an official package. Hence, you can simply do

```julia
add Algencan
```
in package mode in the REPL.

This will download Algencan's code, compile it and make it available to the
Algencan.jl package. **However there is a major caveat here. At this point I am
compiling Algencan without any HSL support. This has major influence on Algencan
behavior and performance. You should use HSL whenever you have access to it.**

You can try to compile Algencan with HSL support you need to get the code from
the [Tango project website](https://www.ime.usp.br/~egbirgin/tango/codes.php)
and compile it yourself, following the authors instructions to use the HSL
libraries and the directions below. Note that the Algencan library has to be
available before installing Algencan.jl so that it can be used by it at
installation. So this compilation needs to be done before doing the
`Pkg.clone`.

### Hints to self compilation of Algencan with HSL libraries

1. Add the option `-fPIC` to  `CFLAGS` and `FFLAGS` in the top of the main
Makefile. Change any numbered compiler version to use the default one in your
system. For example `gcc-4.9` should become `gcc`.

1. Prepare your HSL code as instructed in the `README` file you got from
Algencan. It should be located in `sources\hsl`.

1. Go back to the initial Algencan dir.

1. Type `make` and compile Algencan.

1. Move to the `lib` directory, where you can find the `libalgencan.a` file
and type:
    ```bash
    gcc -shared -o libalgencan.so -Wl,--whole-archive libalgencan.a \\
    -Wl,--no-whole-archive -lgfortran -L$PWD -lhsl
    ```

1. You should now have a file named `libalgencan.so` in the `lib` directory.

1. Create a environmental library named `ALGENCAN_LIB_DIR` pointing to the
`lib` directory. You can proceed to install Algencan.jl as instructed above.
