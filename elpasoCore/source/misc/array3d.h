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

#ifndef INFAM_ARRAY3D_H
#define INFAM_ARRAY3D_H

#include "petsc.h"

//! @brief array used to store values of the shape functions
//! @author Dirk Clasen
//! @date 29.04.2007
//! The copy constructor will only create flat copies of the original
//! object.
class cArray3d {
 private:
  PetscReal* data;  ///< evaluated shapefunctions
  bool isCopy;      ///< used to check, if data only is a copy of another object
  int nn;           ///< number of nodes
  int nf;           ///< number of functions (N  N,xi  N,eta  N,zeta)
  int ng;           ///< number of Gauss points used

 public:
  cArray3d();
  cArray3d(const cArray3d& other);
  ~cArray3d();

  inline int getNumberOfFunctions(void) const { return nf; }
  inline int getNumberOfNodes(void) const { return nn; }
  inline int getNumberOfGaussPoints(void) const { return ng; }

  //! allocate memory
  void initialize(int deriv, int nnod, int n);

  //! return pointer to first entry of the data field (used to "copy" the
  //! object)
  inline double* getBase(void) const { return data; }

  //! access a single value of the data array (read/write)
  inline double& operator()(int deriv, int k, int n) {
    return data[n * nf * nn + deriv * nn + k];
  }

  //! access a single value of the data array (read only)
  inline double operator()(int deriv, int k, int n) const {
    return data[n * nf * nn + deriv * nn + k];
  }
};

#endif
