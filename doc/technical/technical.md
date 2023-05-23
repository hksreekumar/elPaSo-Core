# Developer guide

Welcome to the Developer Guide! Happy to have you onboard and your contribution to improve elPaSo.

The developer circle of elPaSo include the internal researchers and collaborators of the Institute for Acoustics, students and external volunteers. elPaSo-Core is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

elPaSo strives to adhere to the best programming practices like clean code, sustainable research code and efficient software design. Therefore, we encourage the developers to continue this saga.

## Development timeline of elPaSo

As a motivation, elPaSo has a long legacy since 1996. Become a part of it! Following are the important milestones:
|               |               |
|---------------|---------------|
| 1996         | Founded with FEM-BEM coupling <sup>1</sup> |
| 1997          | BEM development and implementation <sup>1</sup>  [H. Antes (1998)]|
| 1997         | Starting point FEM <sup>1</sup> [S. C. Langer, H. Antes (2003)] |
| 2006         | Starting point SBFEM  <sup>1</sup> [L. Lehmann, S. C. Langer, D. Clasen (2006)] |
| 2009         | Communication tools enabled for parallel computations  <sup>1</sup> [M. Schauer et al. (2011)] |
| 2014         | Modeling of flow-induced sound in porous materials <sup>1</sup> [S. C. Beck, S. C. Langer (2014)] |
| 2016         | Uncertainty quantification and sensitivity analysis routines <sup>2</sup> |
| 2019         | Performance optimizations with hybrid parallelization  <sup>3</sup> [H. K. Sreekumar, C. Blech, S. C. Langer (2021)]|
| 2020         | Adapting the fast HDF5 binary file-system <sup>3</sup> |
| 2020         | SURESOFT framework: Development of CI/CD framework, automated testing <sup>3</sup> [H. K. Sreekumar, C. Blech, S. C. Langer (2023)] |
| 2021         | Starting point MOR  <sup>3</sup> |
| 2021         | Fluid structure coupling for non-conforming meshes  <sup>3</sup> |
| 2022         | Starting point parametric MOR  <sup>3</sup> [H. K. Sreekumar et al. (2022,2023)] |
| 2023         | Efficient solving strategies for aircraft cabin noise assessment  <sup>3</sup> [C. Blech, H. K. Sreekumar, Y. Hüpel, S. C. Langer (2023)] |

<sup>1</sup> TU Braunschweig, Institut für Angewandte Mechanik <br>
<sup>2</sup> TU Braunschweig, Institut für Konstruktionstechnik <br>
<sup>3</sup> TU Braunschweig, Institut für Akustik

## Source code

Major amount of source code that belong to the elPaSo solver itself is coded in `C++`. However, elPaSo supporting tools that performs various tasks that includes the pre-/post-visualizing tool, model order reduction toolkit are coded in `python`.

<img src = "../images/sloc.JPG" width="700px">

## Architecture

elPaSo development is targetted towards linux systems so as to finally perform efficient simulations on HPC platforms. The software architecture is depicted in the figure below:

<img src = "../images/elpaso_arch.JPG" width="900px">

## How to Contribute?

You are welcome to contribute to the project. See [Contributing](./contributing.md) for details.

## Maintainance Guide

See documentation on maintaining different entities of elPaSo in [Maintaining](./sustainability/maintainance/maintaining.md).