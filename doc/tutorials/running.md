# Executing the elPaSo Solver

Input to elPaSo is currently limited to tailor-made HDF5 file format produced by the elPaSo Preprocessing tool. Support for other common file formats is in development. However, you can develop a new parser to support your input file by adapting from the parser interface class in cFemParserInterface.h. Contact us if you require support.

## elPaSo Core as a vibroacoustic solver

If you would like to use elPaSo as a vibroacoustic solver, you may need to choose one of the following executable according to requirement:
1. elpaso: Executable compiled with real data-types and hence computations without complex domain (mainly used for time-domain computations).
2. elpasoC: Executable compiled with complex data type (supports other analysis types like static, frequency-domain computations and many more).

## Execution with elPaSo HDF5 as input

### Simple run

elPaSo can be started by direct use of the compiled files elpaso and elpasoC. The input file is needed in the folder where the calculation is started. A simple serial run with standard settings can be made by the command:

```bash
elpasoC -c -inp myInputFile.hdf5
```

A check of the input file is started by the command: 

```bash
elpasoC -check -inp myInputFile.hdf5
```

To run problems with the executable elpaso having real-only petsc datatype, follow the same commands with the respective executable name:

```bash
elpaso -c -inp myInputFile.hdf5
```

Unit-tests can be performed by running the test executable: 

```bash
elpasoT
```

### Parallel execution

In order to execute elpaso in parallel may some additional MPI commands have to be used. These commands are depending on the computer used, e.g.,

```bash
mpirun -np 4 elpasoC -c -inp myInputFile.hdf5
```

In order to run the program using 4 processes. <tt>mpirun</tt> is the instance supplied from an MPI vendor. For GNU compiler build, we use OpenMPI and for INTEL compiler build, we use Intel MPI.

(hpc-execution)=
### HPC execution 🚀

Prepare the job script and submit to the queing system using <tt>sbatch</tt>:

```bash
sbatch example_job_file.job
```

Some job script templates specific to computation with the TU Braunschweig Phoenix cluster:

1. elPaSo parallel execution with OpenMPI
2. [elPaSo parallel execution with Intel MPI](./job_files/intelmpi.job)
3. elPaSo hybrid MPI+OMP parallel exection with Intel MPI and OMP (make sure that you have the hybrid elPaSo executable)
    - [Jobscript for hybrid MUMPS solver](./job_files/mumpshybrid.job) (<tt>-solver 1</tt>)
    - [Jobscript for hybrid CPARDISO solver](./job_files/cpardisohybrid.job) (<tt>-solver 5</tt>)
4. [Simple sleep job script](./job_files/sleep.job) (use at your responsibility!)

## Execution with ABAQUS INP as input

```{note}
Under construction. Content will be available soon.
```