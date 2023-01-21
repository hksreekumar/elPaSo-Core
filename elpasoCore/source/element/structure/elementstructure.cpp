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

#include "elementstructure.h"

cElementStructure::cElementStructure(
  short NumberOfNodes, short NumberOfDofsPerNode, short NumberOfGaussPoints) :
  cElementFEM(NumberOfNodes, NumberOfDofsPerNode, NumberOfGaussPoints)
{
  // empty
}


cElementStructure::cElementStructure(const cElementStructure &other) :
  cElementFEM(other)
{
  // empty
}


cElementStructure::~cElementStructure()
{
  // empty
}


void cElementStructure::assembleDynamicStiffnessMatrix(const PetscReal &omega, cElementMatrix &EM)
{
  // --- the computation has been modified due to some problems
  //     within poroelastic elements. Its now save to set values
  //     to the element matrices instead of adding them.
  //     The old procedure was:
  //       assembleMassMatrix(EM);
  //       infam::scale(EM, -omega*omega;
  //       assembleStiffnessMatrix(EM);
  cElementMatrix lK( EM.rows(), EM.cols() );
  cElementMatrix lM( EM.rows(), EM.cols() );

  assembleMassMatrix(lM);
  infam::scale(lM, -omega*omega);
  assembleStiffnessMatrix(lK, NULL , NULL);

  add(lM, EM);
  add(lK, EM);
}


void cElementStructure::assembleDynamicLoadVector(const PetscReal &omega, cElementVector &LV, cElementMatrix &EM)
{
  assembleLoadVector( LV, EM, NULL, NULL );
}
