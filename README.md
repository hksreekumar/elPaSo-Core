<img src="./doc/logo/InA-TECH-elPaSo-rgb.png" align="right" width=100 height=auto/>

# elPaSo Core - Elementary Parallel Solver

Readme file contains necessary information for using and developing the code project. The contents are structured as per the following sections:

		[SECTION 01 - Introduction]
		[SECTION 02 - Documentations]
		[SECTION 03 - Usage]
		[SECTION 04 - Software dependencies]
		[SECTION 05 - Testing]
		[SECTION 06 - Downloads & Updates]
		[SECTION 07 - License]
		[SECTION 08 - Contributing]
		[SECTION 09 - Support]
		
## Introduction

Many fields of engineering are related to vibration analysis, modal analysis and wave propagation (including unbounded propagation).The elPaSo parallel solver is suitable for acoustic and structural analysis and has been under constant development for many years. The code is based on FEM, BEM, and SBFEM and comprises various material models and element types.

The code is permanently improved. Cooperation with new partners for further challenges is welcome.

## Documentations
Various documentations are available for using/developing elPaSo.
     Doxygen report                  : ./doc/doxygen
     Manual                          : ./doc/manual

## Usage

elPaSo Core can be used as a standalone software and can also be used as a FEM library in other code projects. The software executable and the library can be used by building and compiling elPaSo locally (follow installation instruction in the next section) or can be used with the of docker container (information coming soon).

### Running executable in docker

     # run docker
     docker run -it git.rz.tu-bs.de:4567/akustik/elpaso-core/elpasocore:latest
     # within docker run
     elpasoC -c -inp filename.hdf5
 
## Software dependencies

	 Software Requirements:
        Coreform Cubit          (Stable: 2021.4) - Modelling tool for FEM model preparation
        elPaSo Pre-/Post-Tool         - Tool to produce elPaSo input files
        Microsoft Visual Studio Code  - IDE for code development (with docker containers)
        Cmake                   (Stable: 3.16.3) - Project building
        GNU/Intel compilers     (Stable: GNU 8)/(Intel parallel studio 2020)
        PETSc                   (Stable: 3.16.1)
        SLEPc                   (Stable: 3.16.0)
        ARPACK                  (Stable: 2.1)
        Intel MKL               (Stable: Intel parallel studio 2020)
        HDF5                    (Stable: 1.12.0)
        OpenMPI/IntelMPI        (Stable: OpenMPI 3.1.6/Intel parallel studio 2020)
        
 
## Testing

elPaSo Core uses GoogleTest for running its unit test.
### Running tests in docker
     # run docker
     docker run -it git.rz.tu-bs.de:4567/akustik/elpaso-core/elpasocore:latest
     # within docker run
     elpasoT

## Downloads & Updates

The most recent version of elPaSo can be found at the following locations: (a) Repository: https://git.rz.tu-bs.de/akustik/elPaSo-Core

## License

elPaSo-Core is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 elPaSo-Core is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with elPaSo-Core (COPYING.txt and COPYING.LESSER.txt). If not, see <https://www.gnu.org/licenses/>.

## Contributing

See CONTRIBUTING.md.

## Support

For support please contact us via email at ina(at)tu-braunschweig.de.
	


