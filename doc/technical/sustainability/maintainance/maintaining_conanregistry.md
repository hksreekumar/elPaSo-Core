# Package Registry
{bdg-secondary}`stable | Conan version 1.55.0`

With the GitLab's package registry feature, a pre-built version of various elPaSo dependencies can be hosted. As a result, elPaSo installation effort is significantly reduced thereby avoiding the error-prone installation of its dependencies. We use the [CONAN](https://conan.io/) package manager to handle our C++ external libraries.

If you would like to maintain these packages, you will need an intial setup of conan in your local system and in addition some package specific information which are described below.

(initial-setup)=
## CONAN setup
Refer to [CONAN documentation](https://docs.conan.io) for basic understanding. Follow the steps below for the intial setup:
- Local installation with python (check for stable conan versions):
    ```bash
    python3 -m pip install conan==1.55.0 --user
    conan --version # if fails, open a new terminal and try again
    ```
- So as for CONAN to identify certain compiler version, use the following for updating your conan settings.yml file:
    ```bash
    wget --no-check-certificate https://cloudstorage.tu-braunschweig.de/dl/fi4KjSQYaym7wdoWEe7eyKi6/conan-settings.yml
    mv conan-settings.yml ~/.conan/settings.yml
    ```
- Setup remote repositories
    ```bash
    conan remote add gitlab https://git.rz.tu-bs.de/api/v4/projects/2591/packages/conan
    conan remote remove conancenter
    conan remote remove conan-center
    ```
 - (not required for elpaso) In case of private repositories, you need to enable access to the registry by generating a token (repository > repostory/user settings > Access Tokens > Generate one with minimum role `developer` with scope `read_api` or if you desire to push packages you require scope `api`):
    ```bash
    conan user gitlab-core-user -r gitlab -p <TOKEN>
    ```

### CONAN Cheatsheet

```bash
conan --version
conan remote list
conan user list
conan search # lists existing packages

# example for installation
conan install petsc-complex/3.16.1@ina+elpaso/stable --build missing --settings build_type=Release --settings compiler=gcc --settings compiler.version=8 --settings compiler.libcxx=libstdc++11
```

## Maintaining packages

Maintaining packages would require the above setup and rights to push packages into the repository. Follow [Initial setup](initial-setup) for basic requirements.

- For the desired package, perform a local installation according to steps in [Package Installation Instructions](./maintaining_packages.md).
- We export a conan package out of the local installation. The individual commands to perform the packaging and uploading to the Gitlab's package registry can be found below according to the respective elPaSo dependencies.

```{warning}
Make sure your system environment is compatible with the elPaSo architechture. Recheck for compiler (GNU or Intel) and package versions.
```

### PETSc-Real

{bdg-secondary}`stable | PETSc 3.16.1`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export PETSC_VERSION=3.16.1

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/petsc-3.16.1
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/conan-petsc/conanfile-petsc$PETSC_VERSION-real-gnu.py petsc-real/$PETSC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f

# Intel compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsI20-PS16/petsc-3.16.1
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/conan-petsc/conanfile-petsc$PETSC_VERSION-real-intel.py petsc-real/$PETSC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload petsc-real/$PETSC_VERSION@ina+elpaso/stable --all -r gitlab
```

### PETSc-Complex

{bdg-secondary}`stable | PETSc 3.16.1`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export PETSC_VERSION=3.16.1

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/petsc-3.16.1
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/conan-petsc/conanfile-petsc$PETSC_VERSION-complex-gnu.py petsc-complex/$PETSC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f

# Intel compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsI20-PS16/petsc-3.16.1
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/conan-petsc/conanfile-petsc$PETSC_VERSION-complex-intel.py petsc-complex/$PETSC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload petsc-complex/$PETSC_VERSION@ina+elpaso/stable --all -r gitlab
```

### SLEPc-Real

{bdg-secondary}`stable | SLEPc 3.16.0`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export SLEPC_VERSION=3.16.0

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/slepc-3.16.0
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-slepc/conanfile-slepc$SLEPC_VERSION-real-gnu.py slepc-real/$SLEPC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f

# Intel compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsI20-PS16/slepc-3.16.0
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-slepc/conanfile-slepc$SLEPC_VERSION-real-intel.py slepc-real/$SLEPC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload slepc-real/$SLEPC_VERSION@ina+elpaso/stable --all -r gitlab
```

### SLEPc-Complex

{bdg-secondary}`stable | SLEPc 3.16.0`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export SLEPC_VERSION=3.16.0

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/slepc-3.16.0
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-slepc/conanfile-slepc$SLEPC_VERSION-complex-gnu.py slepc-complex/$SLEPC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f

# Intel compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsI20-PS16/slepc-3.16.0
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-slepc/conanfile-slepc$SLEPC_VERSION-complex-intel.py slepc-complex/$SLEPC_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload slepc-complex/$SLEPC_VERSION@ina+elpaso/stable --all -r gitlab
```

### ARPACK-Real

{bdg-secondary}`stable | ARPACK 2.1`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export ARPACK_VERSION=2.1

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/ARPACK
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-arpack/conanfile-arpack$ARPACK_VERSION-real-gnu.py arpack-real/$ARPACK_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f

# Intel compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsI20-PS16/ARPACK
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-arpack/conanfile-arpack$ARPACK_VERSION-real-intel.py arpack-real/$ARPACK_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload arpack-real/$ARPACK_VERSION@ina+elpaso/stable --all -r gitlab
```

### ARPACK-Complex

{bdg-secondary}`stable | ARPACK 2.1`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export ARPACK_VERSION=2.1

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/ARPACK
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-arpack/conanfile-arpack$ARPACK_VERSION-complex-gnu.py arpack-complex/$ARPACK_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f

# Intel compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsI20-PS16/ARPACK
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-arpack/conanfile-arpack$ARPACK_VERSION-complex-intel.py arpack-complex/$ARPACK_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload arpack-complex/$ARPACK_VERSION@ina+elpaso/stable --all -r gitlab
```

### OpenMPI

{bdg-secondary}`stable | OpenMPI 3.1.6`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export OPENMPI_VERSION=3.1.6

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/openmpi-3.1.6
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-openmpi/conanfile-openmpi$OPENMPI_VERSION-gnu.py openmpi/$OPENMPI_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f
conan upload openmpi/$OPENMPI_VERSION@ina+elpaso/stable --all -r gitlab
```

### HDF5

{bdg-secondary}`stable | HDF5 1.12.0`

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_research_module_repository>/elpasoResearch/envci/conanPackages
export HDF5_VERSION=1.12.0

# GNU compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsGNU-PS16/hdf5-1.12.0
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoResearch/envci/conanPackages/conan-hdf5/conanfile-hdf5$HDF5_VERSION-gnu.py hdf5/$HDF5_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f

# Intel compiled
export PATH_LOCAL_LIB=/home/sreekumar/software/libsconan/PS16/libsI20-PS16/hdf5-1.12.0
cd $LOCAL_LIB_PATH
conan export-pkg $PATH_CONAN_ELPASO_REPOS/NEW_elpaso-Research/elpasoResearch/envci/conanPackages/conan-hdf5/conanfile-hdf5$HDF5_VERSION-intel.py hdf5/$HDF5_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload hdf5/$HDF5_VERSION@ina+elpaso/stable --all -r gitlab
```

### elPaSo Core

{bdg-secondary}`stable | elPaSo Core latest`

#### Executable

```bash
cd <REPOSITORY>/bin
conan export-pkg ../elpasoCore/conan/conan-elpasocore/conanfile-elpasocore-gnu.py elpasocore/23.05.1@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f
conan upload elpasocore/23.05.1@ina+elpaso/stable --all -r gitlab
```

#### Library

```bash
export PATH_CONAN_ELPASO_REPOS=<elpaso_core_module_repository>
export ELPASOCORE_VERSION=23.01.1


# Real
cd $PATH_CONAN_ELPASO_REPOS/bin
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoCore/conan/conan-elpasocore/conanfile-elpasocore-real-gnu.py elpasocore-real/$ELPASOCORE_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoCore/conan/conan-elpasocore/conanfile-elpasocore-complex-gnu.py elpasocore-complex/$ELPASOCORE_VERSION@ina+elpaso/stable -s os=Linux -s compiler=gcc -s compiler.version=8 -s compiler.libcxx=libstdc++11 -f
conan upload elpasocore-complex/$ELPASOCORE_VERSION@ina+elpaso/stable --all -r gitlab

# Complex
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoCore/conan/conan-elpasocore/conanfile-elpasocore-real-intel.py elpasocore-real/$ELPASOCORE_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan export-pkg $PATH_CONAN_ELPASO_REPOS/elpasoCore/conan/conan-elpasocore/conanfile-elpasocore-complex-intel.py elpasocore-complex/$ELPASOCORE_VERSION@ina+elpaso/stable -s os=Linux -s compiler=intel -s compiler.version=19.1 -s compiler.libcxx=libstdc++11 -f
conan upload elpasocore-real/$ELPASOCORE_VERSION@ina+elpaso/stable --all -r gitlab
```