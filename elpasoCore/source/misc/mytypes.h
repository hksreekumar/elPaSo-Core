/* Copyright (c) 2023. Authors listed in AUTHORS.md

 * This file is part of elPaSo-Core.

 * elPaSo-Core is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your option)
 * any later version.

 * elPaSo-Core is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.

 * You should have received a copy of the GNU Lesser General Public License along
 * with elPaSo-Core (COPYING.txt and COPYING.LESSER.txt). If not, see
 * <https://www.gnu.org/licenses/>. 
 */

//! @brief definition of global used datatypes
//! @author Dirk Clasen
//! @date 26.09.2005
//! Within this file all globally uses datatypes are defined.
//! Also most of the systemheaders are included here.
//! petscconf.h defines all stuff that belongs to PETSC like e.g.
//! complex numbers (PETSC_USE_COMPLEX).

#ifndef INFAM_MYTYPES_H
#define INFAM_MYTYPES_H

#include <algorithm>
#include <bitset>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <string>
#include <valarray>
#include <vector>

#include "../basics/exceptions/exceptions.h"
#include "../math/mathlibrary.h"
#include "../math/matrix.h"
#include "../math/vektor.h"
#include "./log/logging.h"

//#include "../../../basics/misc.h"

//! this is the suffix of the input file.
//! If you change it here, it will be changed for
//! the whole code.
#if (defined HAVE_HDF5 && defined USE_HDF5_INP)
const std::string cstInputFileSuffix = ".hdf5";
#else
const std::string cstInputFileSuffix = ".ak3";
#endif  // HAVE_HDF5

//! this enumeration is used by the elements to provide information
//! of their shape. It is evaluated by the postprocessor when he
//! wants to export the mesh or the results to a Tecplot DAT file.
//! Tecplot needs for every element type a different ZONE. In
//! combination with ePhysicsType the elements can be categorized
//! adequately.
enum eElementShape {
  Point,
  Beam,
  Tria,
  Quadrilateral,
  QuadSerendipity,
  Tetrahedron,
  Hexahedron
};

const short cstNumberOfStressDofs = 9;

//! used for adressing stresses and
//! also the components of grad(p) (for 3d poroelastic elements)
enum eStresses {
  sigma11,
  sigma22,
  sigma33,
  sigma23,
  sigma13,
  sigma12,
  gradpx,
  gradpy,
  gradpz
};

//! this enumeration will be used by the elements to specify their element
//! type. This is much easier to use than nested dynamic_cast<> calls.
enum eElementType {
  Unknown,
  Pointmass,
  Spring,
  SpringBC,
  Cable,
  Continuum,
  TimoshenkoBeam,
  Membrane,
  MembraneDrilling,
  MindlinPlate,
  KirchhoffPlate,
  KienzlerPlate,
  PlaneShell,
  PlaneShellDrilling,
  Helmholtz,
  PoroTrash,
  Poro3dUP,
  BeamAcoustic,     // coupling element beam and acoustic fluid
  ShellAcoustic,    // coupling element plane shell or plate and acoustic fluid
  MindlinPoro3dUP,  // coupling element Mindlin plate / 3d poroel. Element u-p
  PoroShellAcoustic,  // coupling element 2d poro elements and acoustic fluid
  MindlinEquiporo,    // coupling element Mindlin plate / equiporo fluid
  SwebemBnd,
  FFPoro3dUP,         // coupling element fluid flow / 3d poroel. Element u-p
  FFShell,            // coupling element fluid flow with shell
  PlaneShellPoro2dUP  // coupling element poro 2d (poro shell Kienzler) /
                      // elastic plane shell
};

//! this information will be given by the elements and is used
//! to split the domains when the results are exported to
//! Tecplot DAT files
enum ePhysicsType {
  Undefined = 0,
  Acoustics = 1,
  Elastodynamics = 2,
  FemBem = 3
};

//! cMatrix and cVector are datatypes for real valued tensors
//! e.g. the Jacobian. The stiffnessmatrix and the massmatrix
//! can consist of complex numbers if we use the frequency domain
//! analysis. Therefore, we have to specify two different types
//! for them: cElementMatrix and cElementVector
typedef infam::tMatrix<PetscReal, PetscInt> cMatrix;
typedef infam::tVector<PetscReal, PetscInt> cVector;

typedef infam::tMatrix<PetscScalar, PetscInt> cElementMatrix;
typedef infam::tVector<PetscScalar, PetscInt> cElementVector;

//! geometrical epsilon
const PetscReal cstGeomEps = 1.0e-5;

/**
 * mathematical constants
 */
const PetscReal eulerNumber = 2.71828182845904523536;  //< euler number
const PetscReal pi = 3.14159265358979323846;           ///< \f$ \pi \f$
const PetscReal one_over_three =
    0.33333333333333333333;                             ///< \f$ \frac{1}{3} \f$
const PetscReal one_over_six = 0.16666666666666666667;  ///< \f$ \frac{1}{6} \f$
const PetscReal one_over_nine =
    0.11111111111111111111;                            ///< \f$ \frac{1}{9} \f$
const PetscReal one_over_27 = 0.03703703703703703703;  ///< \f$ \frac{1}{27} \f$
const PetscReal two_over_three =
    0.66666666666666666667;  ///< \f$ \frac{2}{3} \f$
const PetscReal two_over_nine =
    0.22222222222222222222;  ///< \f$ \frac{2}{9} \f$
const PetscReal four_over_three =
    1.33333333333333333333;  ///< \f$ \frac{4}{3} \f$

const PetscScalar plus_one = 1.0;
const PetscScalar minus_one = -1.0;

// --- the type MPI_UB is not defined when PETSc
//     is compiled in uniprocessor mode (we compile
//     PETSc in this mode only on Windows)
#ifdef _petsc_mpi_uni
#ifndef MPI_UB
#define MPI_UB ((MPI_Datatype)0x4c000011)
#endif
#endif

//! @brief used to leave program on error
//! @author Dirk Clasen
//! @date 17.02.2006
//! if we have a MPI program we need to tell all other
//! processes that an exception occurred - otherwise
//! we'll simply use exit().
#ifdef PETSC_HAVE_MPI
inline void ExitApp(void) {
  cLogging::closeLogFile();
  PetscPrintf(PETSC_COMM_SELF, "Terminating - ExitApp() was called ...\n");
  MPI_Abort(PETSC_COMM_WORLD, PetscGlobalRank);
}
#else
inline void ExitApp(void) {
  cLogging::closeLogFile();
  PetscPrintf(PETSC_COMM_SELF, "Terminating - ExitApp() was called ...\n");
  exit(-1);
}
#endif

template <class T>
inline T infamMaxAbs(const T &a, const T &b) {
  if (std::abs(a) > std::abs(b))
    return a;
  else
    return b;
}

//! check if values of vector vec are a subset of set all
inline bool vectorIsSubset(std::vector<PetscInt> &subset,
                           std::set<PetscInt> &all) {
  std::set<PetscInt>::iterator itset;

  std::sort(subset.begin(), subset.end());

  itset = all.begin();
  for (int k = 0; k < (int)subset.size(); k++) {
    itset = std::find(itset, all.end(), subset[k]);
    if (itset == all.end()) return false;
  }

  return true;
}

// function to sort the vector in abs-order
inline bool abs_order(PetscReal i, PetscReal j) {
  return (std::abs(i) < std::abs(j));
}
#endif
