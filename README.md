# Automaton
This project using 2D decomposition methods to simulation of Game of lifes.

## Compile
Using mpicc to compile the program. To compile on Cirrus, first we should

```
module load mpt 
module load intel-compilers-19
```

Then

```
make
```
The executable file will be named as 'automaton'


## Run

To run on login node

```
mpicc -n <NPROC> ./automaton <SEED> <L> <rho>
```
where NPROC must and SEED must be specified. L and rho are set to default values as 768 and 0.49 if missing to specify. 


## Submission on Cirrus


```
sbatch cirrus_auto_2d.job
```
`<NPROC>` in *.job is specified by ntask.


## Unit Test
To run the unit test, we should compile the unit test by

`make -f UnitTest`

Then using 4 process to run

`mpirun -n 4 ./unit_test`

##  