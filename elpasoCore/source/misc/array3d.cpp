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

#include "array3d.h"

cArray3d::cArray3d()
{
  data   = 0;
  isCopy = false;
  nf = nn = ng = 0;
}


cArray3d::cArray3d(const cArray3d &other) :
  nn( other.getNumberOfNodes() ),
  nf( other.getNumberOfFunctions() ),
  ng( other.getNumberOfGaussPoints() )
{
  // --- just create a flat copy
  //     we know, what we're doing
  data = other.getBase();
  isCopy = true;
}


cArray3d::~cArray3d()
{
  if (isCopy) {
    data = 0;
  }
  else {
    if (data != 0) {
      delete[] data;
      data = 0;
    }
  }
}


void cArray3d::initialize(int deriv, int nnod, int n)
{
  nf = deriv;
  nn = nnod;
  ng = n;

  data = new double[deriv * nnod * n];
  for (int k=0; k<deriv*nnod*n; k++)
    data[k] = 0.;
}
