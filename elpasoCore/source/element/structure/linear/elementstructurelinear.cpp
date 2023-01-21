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

#include "elementstructurelinear.h"

cElementStructureLinear::cElementStructureLinear(
  short NumberOfNodes, short NumberOfDofsPerNode, short NumberOfGaussPoints) :
  cElementStructure(NumberOfNodes, NumberOfDofsPerNode, NumberOfGaussPoints)
{
  m_Material = NULL;
}


cElementStructureLinear::cElementStructureLinear( const cElementStructureLinear &other) :
  cElementStructure(other)
{
  setMaterial( other.getMaterial( ) );
}


cElementStructureLinear::~cElementStructureLinear()
{
  m_Material = NULL;
}


void cElementStructureLinear::setMaterial(cMaterial *ptrMaterial)
{
  if (ptrMaterial == NULL)
    throw cException("***  Pointer to material is NULL  ***", __FILE__, __LINE__);

  m_Material = dynamic_cast<cMaterialStructure *>(ptrMaterial);

  if (m_Material == NULL) {
    message("Id of the material : %d\n", ptrMaterial->getId());
    throw cException("***  wrong type of material for structural element ***", __FILE__, __LINE__);
  }
}
