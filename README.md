# MPI-Game-of-Life
A simple and rough MPI implementation of Conway's Game of Life.

## Prerequisites
* GNU software (GCC 9, Make)
* MPI (tested with Open MPI 4.1.0)
* Python 3.x with numpy and Pillow installed

Developed and tested on a single machine with Ubuntu 20.04 LTS. 

## Compilation
A Makefile has been provided to build and compile the target executable. Simply run `make` in the top level of the project directory.

The target executable `lifeMPI` is created in the same location.

## Usage
On a single machine

```
mpiexec [--use-hwthread-cpus] -np <#> ./lifeMPI <args> [opts]
```
* `--use-hwthreads-cpus` enables use of logical processors (e.g. hyper-threaded cores). This is not likely to improve performance, and may negatively impact performance.

where

```
args :  nrows (uint) --> number of rows in the grid
        ncols (uint) --> number of columns in the grid
        ngens (uint) --> number of generations to simulate
        EITHER filename (string) OR prob_alive (float) --> initial state of the grid

opts :  --periodic --> enable periodic boundaries [default: disabled]
        --dump=<uint> --> num generations between output file dumps [default: 101]
```

Examples:
```
mpiexec -np 4 ./lifeMPI 256 256 500 0.10 --periodic

mpiexec -np 9 ./lifeMPI 247 331 250 grid.bin --dump=63 --periodic
```

**NOTES**
* The options `--periodic` and `--dump=##` must be specified after the 4 input arguments (but the order does not matter).
* The input file must be a raw binary file containing 0's and 1's represented by single bytes.
* The grid dimensions given (`nrows`, `ncols`) must match exactly the input binary file.
* A binary file can be easily generated from a CSV (or any other delimted text file) using a `Python` commnd, e.g.:

```
python -c "import numpy; numpy.loadtxt('grid.csv', dtype=numpy.uint8, delimiter=',').tofile('grid.bin')"
```

## Visualisation

The program will dump binary files of the grid state for every generation to `outfiles/` (directory created during the build process, and can be automatically cleaned with `make clean`). 

An animated gif of the grid state for every generation can be created by running the provided script:
```
python scripts/animate_life.py
```

## Testing

A verified python implemenation of Conway's Game of Life is provided (`scripts/automata.py`) for testing accuracy of the MPI implementation. 

A simple test can be performed by running `make test` (or `python scripts/test_life.py`) which compares the final states of both implementations.
