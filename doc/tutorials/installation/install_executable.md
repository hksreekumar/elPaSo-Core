# Install elPaSo executable

Following are the different ways to install an elPaSo executable. If you would like to get started with elPaSo as easy as possible follow the `recommended` installation routine.

## Ways to install elPaSo executable 🔧

Choose your favourite way to install elPaSo:

1. For **least complicated** and easy to start executable {bdg-info}`easy` {bdg-secondary}`recommended for users`, perform:
    1. [Install a pre-built elPaSo via CONAN package manager](prebuilt-elpaso)
1. For easy installation {bdg-secondary}`recommended for developers/users`, perform:
    1. [Install dependencies via CONAN package manager](dependencies-conan)
    2. [Build elPaSo](build-elpaso)
1. If you have **Docker** support {bdg-secondary}`recommended for users`, you can directly use a latest elPaSo image; perform:
    1. [Install elPaSo Docker](./install_docker.md)
1. For advanced installation {bdg-primary}`recommended for developers/users`, perform:
    1. [Install elPaSo dependencies from scratch](dependencies-scratch)
    2. [Build elPaSo](build-elpaso)

```{important}
Make sure you have a compatible compiler ready. Refer to [External Dependencies](../../technical/sustainability/maintainance/maintaining_libraries.md) for basic OS requirements and compilers. Contact us if you require support compiling with compilers other than GNU.
```

(prebuilt-elpaso)=
## Install a pre-built elPaSo via CONAN package manager

```{note}
Also works in [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install).
```

Make sure you have the CONAN setup according to [CONAN Setup](initial-setup). Then install elPaSo with:

```bash
# requirements
## 1. Compiler gcc/g++/gfortran version >= 8
##        sudo apt-get install gcc g++ gfortran
## 2. Ubuntu version >= 20

# install elPaSo prebuilt
conan install elpasocore/23.05.1@ina+elpaso/stable --build missing --settings build_type=Release --settings compiler=gcc --settings compiler.version=8 --settings compiler.libcxx=libstdc++11 -g virtualrunenv

# prepare environment for correct paths - or add the following in ~/.bashrc
source activate_run.sh
export OPAL_PREFIX=~/.conan/data/elpasocore/23.05.1/ina+elpaso/stable/package/56e0cf6d16ee57367a0661ab743f4e43b29223f8/
export PATH=$PATH:$OPAL_PREFIX
```

```{note}
Currently only GNU build version of elPaSo is available via package manager.
```

### Running in Windows Subsystem for Linux

1. Open a terminal at your desired location in windows
2. Open the ubuntu distribution with command: `wsl`
3. In the wsl ubuntu terminal, set the `~/.bashrc` file with the environment flags mentioned before.
4. Run elpaso by calling `elpasoC -c -inp example.hdf5`.

## Collect dependencies

The section describes how to prepare the required 3rd party dependencies required prior to elPaSo building.

(dependencies-conan)=
## Install dependencies via CONAN package manager

You only need to have the correct CONAN setup according to [CONAN Setup](initial-setup).

(dependencies-scratch)=
## Install elPaSo dependencies from scratch

Follow the various procedures in [External Dependencies](../../technical/sustainability/maintainance/maintaining_libraries.md) for installing the required elPaSo dependencies.

(build-elpaso)=
## Build elPaSo

At this point, you have all the elPaSo requirements prepared. Follow the steps below to install an elPaSo executable:

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
    This creates a config.cmake file with your hostname in the "cmake" folder. The file contains the default configuration. Check for correct paths and versions. If you get the [dependencies via CONAN](dependencies-conan) set <tt>OPTION(LINK_CONAN		"CONAN LIBS"		ON)</tt>, or else if the [dependencies are installed from scratch](dependencies-scratch) set <tt>OPTION(LINK_CONAN		"CONAN LIBS"		OFF)</tt>.
3. Build elPaSo executable
    - With GNU compiler
        ```bash
        cd <REPOS_PATH>
        mkdir project
        make elpaso-gnu
        make elpaso-build
        ```
        When sucessfull, `bin` folder contains the executables.
    - With Intel compiler - Install intel-parallel-studio and configure the ./cmake/\$HOSTNAME.config.cmake with the respective directory where the installation is done with INTEL_DIR and the respective version number to INTELMPI_VERSION.
        ```bash
        cd <REPOS_PATH>
        mkdir project
        make elpaso-intel
        make elpaso-build
        ```
        When sucessfull, `bin` folder contains the executables.
4. You are ready to run the executable. Refer to [Executing the elPaSo Solver](../running.md) for more.
