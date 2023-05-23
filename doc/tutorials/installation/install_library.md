# Install elPaSo as library

```{important}
Make sure you have a compatible compiler ready. Refer to [External Dependencies](../../technical/sustainability/maintainance/maintaining_libraries.md) for basic OS requirements and compilers. Contact us if you require support compiling with compilers other than GNU.
```

## Installation

Follow the steps below:
1. Checkout the current version of the elPaSo Core project
    ```bash
    cd <YOUR_BUILD_PATH>
    git clone https://git.rz.tu-bs.de/akustik/elPaSo-Core
    cd elPaSo-Core
    ```
2. Configure a cmake configuration file for your system.
    ```bash
    cd <REPOS_PATH>/cmake
    cp YOURCOMPUTER.config.cmake ${HOSTNAME}.config.cmake
    ```
    This creates a config.cmake file with your hostname in the "cmake" folder. The file contains the default configuration. To generate the eCore library, active the library generation option in ./cmake/\$HOSTNAME.config.cmake by modifying: <tt>OPTION(GEN_DLIB "elPaSo DYNAMIC LIB"	ON)</tt>
3. Build elPaSo Core library
    - With GNU compiler
        ```bash
        cd <REPOS_PATH>
        mkdir project
        make elpaso-gnu
        make elpasolib-build
        ```
        When sucessfull, elpasoCore-gnu folder contains the necessary entities for linking eCore.
    - With Intel compiler - Install intel-parallel-studio and configure the ./cmake/\$HOSTNAME.config.cmake with the respective directory where the installation is done with INTEL_DIR and the respective version number to INTELMPI_VERSION.
        ```bash
        cd <REPOS_PATH>
        mkdir project
        make elpaso-intel
        make elpasolib-build
        ```
        When sucessfull, elpasoCore-intel folder contains the necessary entities for linking eCore.

## Linking library

Follow the steps below, to link eCore library with CMake. In your FindELPASOCORE.cmake:

1. Include the path <tt>./elpasoCore-<COMPILER>/include</tt> to your cmake project for link to header files
2. Link the eCore dynamic library:
    - For real-valued build, link <tt>./elpasoCore-<COMPILER>/lib/libelpasoCore-intel-cxx-o.so</tt>
    - For complex-valued build, link <tt>./elpasoCore-<COMPILER>/lib/libelpasoCore-intel-cxx-complex.so</tt>
