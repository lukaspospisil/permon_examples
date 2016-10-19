# permon_examples

examples of using PermonQP library (https://github.com/haplav/permon)

Instalation:
- install PETSc library (and set `PETSC_DIR`, `PETSC_ARC`H)
- install PermonQP library (and set `PERMON_DIR`)
- clone repository with examples
```
git clone https://github.com/lukaspospisil/permon_examples
```
- create `build` directory (this directory is included in `.gitignore`)
```
mkdir build
cd build
```
- call cmake to prepare `Makefile`
```
cmake ..
```
- decide which test to build (see the list provided by previous call)
```
cmake -DTEST_01=ON ..
```
- call `make` to compile the code
- if you want to use CUDA, call the cmake with `-DUSE_CUDA=ON`

Usage:
- run compiled example as sequential program
```
./test_01
```
- or run compiled example using standard MPI command
```
mpiexec -n 3 ./test_01
```
