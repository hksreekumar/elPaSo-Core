# External dependencies

elPaSo depend on various 3rd party libraries for state-of-the-art implementations. As already stated, the libraries are compiled with GNU and Intel compilers. Therefore the applicable dependencies are listed below for the two compiler cases:

::::{grid}
:gutter: 2

:::{grid-item-card}
**elPaSo - GNU Build**<br/>
{bdg-secondary}`open source`<br/>
Compiler | GNU compiler version 8<br/>
Dependencies | PETSc, SLEPc, ARPACK, HDF5, OpenMPI
:::

:::{grid-item-card} 
**elPaSo - Intel Build** <br/>
{bdg-info}`recommended` {bdg-info}`performance`<br/>
Compiler | ICC/ICPC version 19.1/Intel Parallel Studio 2020<br/>
Dependencies | PETSc, SLEPc, ARPACK, HDF5, IntelMPI
:::
::::

```{warning}
For intel compiled version, you require a licensed installation of Intel Parallel Studio. Hence, distribution of elPaSo packages compiled with intel is limited. Contact your IT administrator for existing license subscriptions.
```

## Stable setting

|   Dependency                      |   Version                           |
|-------------------------|------------------------------|
|    OS  |  Ubuntu 20.04.5 LTS |
|    GNU/Intel compilers  |   (Stable: GNU 8)/(Intel parallel studio 2020) |
|    cmake                |   (Stable: 3.14.0, Tested: 3.19.4) |
|    PETSc                |   (Stable: 3.19.1) |
|    SLEPc                |   (Stable: 3.19.1) |
|    ARPACK               |   (Stable: 2.1) |
|    Intel MKL            |   (Stable: Intel parallel studio 2020) |
|    HDF5                 |   (Stable: 1.12.0) |
|    OpenMPI/IntelMPI     |   (Stable: OpenMPI 3.1.6/Intel parallel studio 2020) |


## Workstation installation

### Basic requirements
```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends build-essential wget make python libx11-dev zlib1g-dev openssh-server cmake

# versions and path
export CPATH=/opt/elPaSo # installation directory for elpaso
rm -rf $CPATH
mkdir $CPATH


##### COMMON VARIABLES
# PETSc
export PETSC_V=3.19.1
export PETSC_DIR=${CPATH}/petsc-${PETSC_V}
# SLEPc
export SLEPC_V=3.19.1

# HDF5
export HDF5_V=1.12.0


##### FOR GNU BUILD
# GNU Compiler
export GNUC_V=8
# OpenMPI
export OMPI_V=3.1.6
export OMPI_DIR=${CPATH}/openmpi-${OMPI_V}

##### FOR INTEL BUILD
# Intel Compiler
export INTEL_DIR=/software  # installation directory; require sudo rights to install
export INTEL_COMP_V=2020.0.166
export INTEL_MPI_V=2019.6.166
```


### GNU build

1. Compiler
    ```bash
    sudo apt-get install -y --no-install-recommends gcc-${GNUC_V} g++-${GNUC_V} gfortran-${GNUC_V}
    # symbolic linking GNU Compilers
    sudo rm -f /usr/bin/gcc /usr/bin/g++ /usr/bin/gfortran 
    sudo ln -s /usr/bin/gcc-${GNUC_V} /usr/bin/gcc
    sudo ln -s /usr/bin/g++-${GNUC_V} /usr/bin/g++
    sudo ln -s /usr/bin/gfortran-${GNUC_V} /usr/bin/gfortran
    # Check GNU Compiler Versions
    gcc --version
    g++ --version
    gfortran --version
    ```
1. Install OpenMPI
    ```bash
    cd ${CPATH}
    wget --no-check-certificate https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.6.tar.gz 
    tar -xvf openmpi-${OMPI_V}.tar.gz
    rm openmpi-${OMPI_V}.tar.gz
    cd openmpi-${OMPI_V}
    ./configure --prefix=${OMPI_DIR}/gnu-opt
    make -j20 all && make install
    ```
1. Install PETSc
    ```bash
    cd ${CPATH}
    wget --no-check-certificate ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_V}.tar.gz
    tar -zxvf petsc-${PETSC_V}.tar.gz 
    rm petsc-${PETSC_V}.tar.gz
    ## Environmental Variables for PETSc Install
    export PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/bin:${PATH}"
    export LD_LIBRARY_PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/lib" 
    ## Configure and Make routine for PETSc
    cd petsc-${PETSC_V} 
    #
    #
    ###### Install PETSc Complex
    export PETSC_ARCH=gnu-cxx-complex-o 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    ./config/configure.py --PETSC_ARCH=$PETSC_ARCH --with-gnu-compilers=1 --with-shared=1 --with-mpi=1 --with-mpi-dir=${OMPI_DIR}/gnu-opt --PETSC_DIR=${PETSC_DIR} --download-elemental --download-plapack --download-blacs --download-fblaslapack --download-scalapack --download-mumps --download-parmetis --download-metis --download-metis-cmake-arguments=-DMETIS_USE_DOUBLEPRECISION=1 --with-clanguage=cxx --with-scalar-type=complex --with-debugging=0 --with-cxx-dialect=C++11 
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH -j20 all 
    #make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH test 
    #
    #
    ###### Install PETSc Real
    export PETSC_ARCH=gnu-cxx-o 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    ./config/configure.py --PETSC_ARCH=$PETSC_ARCH --with-gnu-compilers=1 --with-shared=1 --with-mpi=1 --with-mpi-dir=${OMPI_DIR}/gnu-opt --PETSC_DIR=${PETSC_DIR} --download-elemental --download-plapack --download-blacs --download-fblaslapack --download-scalapack --download-mumps --download-parmetis --download-metis --download-metis-cmake-arguments=-DMETIS_USE_DOUBLEPRECISION=1 --with-clanguage=cxx --with-scalar-type=real --with-debugging=0 --with-cxx-dialect=C++11 
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH -j20 all 
    #make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH test     
    ```
1. Install ARPACK and SLEPc - complex setting
    ```bash
    cd ${CPATH} 
    #
    #
    #### ARAPACK Installation
    ## Environmental Variables for PETSc Install
    export PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/bin:${PATH}" 
    export LD_LIBRARY_PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/lib" 
    export PETSC_ARCH=gnu-cxx-complex-o 
    export ARPACK_DIR="${CPATH}/ARPACK" 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    # Download ARPACK from repository
    wget --no-check-certificate https://cloudstorage.tu-braunschweig.de/dl/fiGkUYwu3s4DS4NJPxXx6BcC/arpack_c.tar.gz 
    tar -xvf arpack_c.tar.gz 
    rm arpack_c.tar.gz 
    # Basic environmental flags
    export FC=$OMPI_DIR/gnu-opt/bin/mpif77 
    export FC_FLAGS=-fPIC 
    echo $FC && echo $FC_FLAGS 
    # compile ARPACK with new gcc/g++
    cd $ARPACK_DIR 
    make cleanlib 
    make all 
    mkdir $PETSC_ARCH &&
    mkdir $PETSC_ARCH/lib &&
    mv libarpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libarpack_SUN4.a &&
    mv libparpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libparpack_MPI-SUN4.a &&
    #
    #
    #### SLEPc Installation
    cd ${CPATH} &&
    ## Environmental Variables for PETSc Install
    export PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/bin:${PATH}" &&
    export LD_LIBRARY_PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/lib" &&
    export PETSC_ARCH=gnu-cxx-complex-o &&
    export ARPACK_DIR="${CPATH}/ARPACK" &&
    export SLEPC_DIR="${CPATH}/slepc-${SLEPC_V}" &&
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH &&
    ## Download
    wget --no-check-certificate https://slepc.upv.es/download/distrib/slepc-${SLEPC_V}.tar.gz &&
    tar -xvf slepc-${SLEPC_V}.tar.gz &&
    rm slepc-${SLEPC_V}.tar.gz &&
    ## Install
    cd $SLEPC_DIR &&
    ./config/configure.py --with-arpack --with-arpack-dir=$ARPACK_DIR/$PETSC_ARCH &&
    make SLEPC_DIR=$SLEPC_DIR PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH &&
    #make check
    ```
1. Install ARPACK and SLEPc - real setting
    ```bash
    cd ${CPATH} &&
    #
    #
    #### ARAPACK Installation
    ## Environmental Variables for PETSc Install
    export PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/bin:${PATH}"  &&
    export LD_LIBRARY_PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/lib"  &&
    export PETSC_ARCH=gnu-cxx-o  &&
    export ARPACK_DIR="${CPATH}/ARPACK"  &&
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH  &&
    # Basic environmental flags
    export FC=$OMPI_DIR/gnu-opt/bin/mpif77  &&
    export FC_FLAGS=-fPIC  &&
    echo $FC && echo $FC_FLAGS &&
    # compile ARPACK with new gcc/g++
    cd $ARPACK_DIR  &&
    make cleanlib  &&
    make all  &&
    mkdir $PETSC_ARCH &&
    mkdir $PETSC_ARCH/lib &&
    mv libarpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libarpack_SUN4.a &&
    mv libparpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libparpack_MPI-SUN4.a &&
    #
    #
    #### SLEPc Installation
    cd ${CPATH} && 
    ## Environmental Variables for PETSc Install
    export PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/bin:${PATH}" && 
    export LD_LIBRARY_PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/lib" && 
    export PETSC_ARCH=gnu-cxx-o  &&
    export ARPACK_DIR="${CPATH}/ARPACK"  &&
    export SLEPC_DIR="${CPATH}/slepc-${SLEPC_V}"  &&
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH  &&
    ## Install
    cd $SLEPC_DIR  &&
    ./config/configure.py --with-arpack --with-arpack-dir=$ARPACK_DIR/$PETSC_ARCH &&
    make all
    ```
1. HDF5 
    ```bash
    cd ${CPATH} &&
    wget --no-check-certificate https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/CMake-hdf5-${HDF5_V}.tar.gz  &&
    tar -xvf CMake-hdf5-${HDF5_V}.tar.gz  &&
    rm CMake-hdf5-${HDF5_V}.tar.gz  &&
    cd CMake-hdf5-${HDF5_V}/hdf5-${HDF5_V}  &&
    export HDF5_DIR="${CPATH}/hdf5-${HDF5_V}"  &&
    export PATH="${CPATH}/openmpi-${OMPI_V}/gnu-opt/bin:${PATH}"  &&
    ls ${CPATH}/openmpi-${OMPI_V}/gnu-opt/bin  &&
    CC=mpicc CXX=mpicpc ./configure --prefix=${HDF5_DIR}/gnu-opt --enable-parallel --enable-shared --enable-build-mode=production --disable-hl  &&
    make  &&
    make install  &&
    cd ${CPATH}  &&
    rm -r CMake-hdf5-${HDF5_V}
    ```

### Intel build

1. Requirements
    ```bash
    sudo apt-get install -y --no-install-recommends cpio libpango-1.0-0 libasound2 libgtk2.0-dev libfabric-dev
    source ${INTEL_DIR}/intel/compilers_and_libraries_${INTEL_COMP_V}/linux/bin/compilervars.sh -arch intel64 -platform linux
    ```
1. Compiler installation instruction for TU Braunschweig users. For external user, contact IT administrator for valid license.
    ```bash
    cd $CPATH # or any other temporary directory
    wget --no-check-certificate https://campus-software.tu-braunschweig.de/campus-software/Intel/2020/Intel%20parallel%20Studio%20Studio%20XE%202020%20Cluster%20Edition/Linux/parallel_studio_xe_2020_cluster_edition.tgz 
    tar -xvf parallel_studio_xe_2020_cluster_edition.tgz
    rm parallel_studio_xe_2020_cluster_edition.tgz
    #
    cd /opt/parallel_studio_xe_2020_cluster_edition
    # silent configuration file
    wget --no-check-certificate https://cloudstorage.tu-braunschweig.de/dl/fi43EFZAypi3FQGkozURnA2F/elPaSo-Intel2020.cfg
    mv elPaSo-Intel2020.cfg silent.cfg
    ./install.sh --silent silent.cfg
    cd /opt
    rm -r parallel_studio_xe_2020_cluster_edition
    ## Symbolic link intel compilers
    ln -s ${INTEL_DIR}/intel/compilers_and_libraries_${INTEL_COMP_V}/linux/bin/intel64/icc /usr/bin/icc
    ln -s ${INTEL_DIR}/intel/compilers_and_libraries_${INTEL_COMP_V}/linux/bin/intel64/icpc /usr/bin/icpc
    ln -s ${INTEL_DIR}/intel/compilers_and_libraries_${INTEL_COMP_V}/linux/bin/intel64/ifortran /usr/bin/ifortran
    ```
1. Install PETSc
    ```bash
    source ${INTEL_DIR}/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux

    # Download
    cd ${CPATH} 
    wget --no-check-certificate ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PETSC_V}.tar.gz 
    tar -zxvf petsc-${PETSC_V}.tar.gz 
    rm petsc-${PETSC_V}.tar.gz 
    #
    #
    ## Install petsc complex
    export PETSC_ARCH=intel-cxx-complex-o 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    ## Configure and Make routine for PETSc
    cd petsc-${PETSC_V} 
    ./config/configure.py --with-cc=mpiicc --with-cxx=mpiicpc --with-fc=mpiifort --PETSC_ARCH=$PETSC_ARCH --PETSC_DIR=$PETSC_DIR --with-blas-lapack-dir=$INTEL_DIR/intel/mkl/lib/intel64 --with-mkl_pardiso-dir=$INTEL_DIR/intel/mkl --with-mkl_cpardiso-dir=$INTEL_DIR/intel/mkl --download-elemental --download-scalapack --download-mumps --download-parmetis --download-metis --download-metis-cmake-arguments=-DMETIS_USE_DOUBLEPRECISION=1 --with-scalar-type=complex --with-cxx-dialect=C++11 --with-debugging=0 
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH -j20 all 
    #make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH test 
    #
    #
    ## Install petsc real
    export PETSC_ARCH=intel-cxx-o 
    ./config/configure.py --with-cc=mpiicc --with-cxx=mpiicpc --with-fc=mpiifort --PETSC_ARCH=$PETSC_ARCH --PETSC_DIR=$PETSC_DIR --with-blas-lapack-dir=$INTEL_DIR/intel/mkl/lib/intel64 --with-mkl_pardiso-dir=$INTEL_DIR/intel/mkl --with-mkl_cpardiso-dir=$INTEL_DIR/intel/mkl --download-elemental --download-scalapack --download-mumps --download-parmetis --download-metis --download-metis-cmake-arguments=-DMETIS_USE_DOUBLEPRECISION=1 --with-scalar-type=real --with-cxx-dialect=C++11 --with-debugging=0 
    make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH -j20 all 
    #make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=$PETSC_ARCH test 
    ```
1. Install ARPACK and SLEPc - complex setting
    ```bash
    source ${INTEL_DIR}/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
    #
    #
    #### ARAPACK Installation
    cd ${CPATH} 
    export PETSC_ARCH=intel-cxx-complex-o 
    export ARPACK_DIR="${CPATH}/ARPACK" 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    # Download ARPACK from repository
    wget --no-check-certificate https://cloudstorage.tu-braunschweig.de/dl/fiGkUYwu3s4DS4NJPxXx6BcC/arpack_c.tar.gz 
    tar -xvf arpack_c.tar.gz 
    rm arpack_c.tar.gz 
    # Basic environmental flags
    export FC=${INTEL_DIR}/intel/impi/${INTEL_MPI_V}/intel64/bin/mpiifort 
    export FC_FLAGS=-fPIC 
    echo $FC && echo $FC_FLAGS 
    # compile ARPACK
    cd $ARPACK_DIR 
    make cleanlib 
    make all 
    mkdir $PETSC_ARCH
    mkdir $PETSC_ARCH/lib
    mv libarpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libarpack_SUN4.a
    mv libparpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libparpack_MPI-SUN4.a
    #
    #
    #### SLEPc Installation
    cd ${CPATH} 
    ## Environmental Variables for SLEPc Install
    export PETSC_ARCH=intel-cxx-complex-o 
    export ARPACK_DIR="${CPATH}/ARPACK" 
    export SLEPC_DIR="${CPATH}/slepc-${SLEPC_V}" 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    ## Download
    wget --no-check-certificate https://slepc.upv.es/download/distrib/slepc-${SLEPC_V}.tar.gz 
    tar -xvf slepc-${SLEPC_V}.tar.gz 
    rm slepc-${SLEPC_V}.tar.gz 
    ## Install
    cd $SLEPC_DIR 
    ./config/configure.py --with-arpack --with-arpack-dir=$ARPACK_DIR/$PETSC_ARCH
    make all 
    make check
    ```
1. Install ARPACK and SLEPc - real setting
    ```bash
    source ${INTEL_DIR}/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
    #
    #
    #### ARAPACK Installation
    cd ${CPATH} 
    ## Environmental Variables for PETSc Install
    export PETSC_ARCH=intel-cxx-o 
    export ARPACK_DIR="${CPATH}/ARPACK" 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    # Download ARPACK from repository
    wget --no-check-certificate https://cloudstorage.tu-braunschweig.de/dl/fiGkUYwu3s4DS4NJPxXx6BcC/arpack_c.tar.gz 
    tar -xvf arpack_c.tar.gz 
    rm arpack_c.tar.gz 
    # Basic environmental flags
    export FC=${INTEL_DIR}/intel/impi/${INTEL_MPI_V}/intel64/bin/mpiifort 
    export FC_FLAGS=-fPIC 
    echo $FC && echo $FC_FLAGS 
    # compile ARPACK
    cd $ARPACK_DIR 
    make cleanlib 
    make all 
    mkdir $PETSC_ARCH
    mkdir $PETSC_ARCH/lib
    mv libarpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libarpack_SUN4.a
    mv libparpack_$PETSC_ARCH.a $PETSC_ARCH/lib/libparpack_MPI-SUN4.a
    #
    #
    #### SLEPc Installation
    cd ${CPATH} 
    ## Environmental Variables for SLEPc Install
    export PETSC_ARCH=intel-cxx-o 
    export ARPACK_DIR="${CPATH}/ARPACK" 
    export SLEPC_DIR="${CPATH}/slepc-${SLEPC_V}" 
    echo $LD_LIBRARY_PATH && echo ${PETSC_DIR} && echo $PETSC_ARCH 
    ## Download
    wget --no-check-certificate https://slepc.upv.es/download/distrib/slepc-${SLEPC_V}.tar.gz 
    tar -xvf slepc-${SLEPC_V}.tar.gz 
    rm slepc-${SLEPC_V}.tar.gz 
    ## Install
    cd $SLEPC_DIR 
    ./config/configure.py --with-arpack --with-arpack-dir=$ARPACK_DIR/$PETSC_ARCH
    make all 
    make check
    ```
1. HDF5 
    ```bash
    source ${INTEL_DIR}/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux
    #
    #
    cd ${CPATH} 
    wget --no-check-certificate https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/CMake-hdf5-${HDF5_V}.tar.gz 
    tar -xvf CMake-hdf5-${HDF5_V}.tar.gz 
    rm CMake-hdf5-${HDF5_V}.tar.gz 
    cd CMake-hdf5-${HDF5_V}/hdf5-${HDF5_V} 
    export HDF5_DIR="${CPATH}/hdf5-${HDF5_V}" 
    CC=mpiicc CXX=mpiicpc ./configure --prefix=${HDF5_DIR}/intel-opt --enable-parallel --enable-shared --enable-build-mode=production --disable-hl 
    make 
    make install 
    cd ${CPATH} 
    rm -r CMake-hdf5-${HDF5_V}
    ```

## Cluster installation

Following specifics are for installation in TU Braunshcweig Phoenix HPC cluster.

### GNU build

Follow the procedure for local installation except for the following:
1. Loading gnu libraries (installation not required!):
```bash
module load comp/gcc/8.3.0 
module load comp/cmake/3.14.0
```

### Intel build

Follow the procedure for local installation except for the following:
1. Loading intel libraries (installation not required!):
```bash
module load intel-studio-2020
source /cluster/share/intelenv.sh
module load comp/cmake/3.14.0
```

## Other tools

### Code quality linter

We use `clang-tidy` for code quality checks. Install script:

```bash
sudo apt-get install -y --no-install-recommends clang-tidy && \
clang-tidy --version
```

### Python

Python installation is not required for normal use. But following are the stable setting for CI pipeline:


```bash
sudo apt-get install -y --no-install-recommends python python3-pip git
sudo python3 -m pip install tikzplotlib fpdf PyPDF2 python-gitlab gitpython pandas
```

## Older records

Installation for PETSc 3.16.1 package - [See file](./older_/maintaining_libraries-3161.md)