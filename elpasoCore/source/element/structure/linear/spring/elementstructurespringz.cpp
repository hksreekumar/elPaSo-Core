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

#include "elementstructurespringz.h"

cElementStructureSpringz::cElementStructureSpringz() :
  cElementStructureLinear(2, 1, 0)
{
  // empty
}


cElementStructureSpringz::cElementStructureSpringz(const cElementStructureSpringz &other) :
  cElementStructureLinear(other)
{
  // empty
}


cElementStructureSpringz::~cElementStructureSpringz()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureSpringz::getDofs(void) const
{
  std::vector<eKnownDofs> res(1);

  res[0] = disp_x3;

  return res;
}


void cElementStructureSpringz::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  const PetscScalar Cz = m_Material->getCz();
  KM(0,0) =  Cz;   KM(0,1) =  -Cz;
  KM(1,0) = -Cz;   KM(1,1) =  Cz;

  // std::cout << KM << std::endl;
}



void cElementStructureSpringz::assembleMassMatrix(cElementMatrix &MM)
{
  // nothing to do here
}


void cElementStructureSpringz::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
   // nothing to do here
}


std::ostream& cElementStructureSpringz::write(std::ostream &os) const
{
  os << "spring element" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k=0; k<getNumberOfNodes(); k++)
    os << *(m_Nodes[k]);

  return os;
}


std::ostream& cElementStructureSpringz::writeXml(std::ostream &os) const
{
  os << "<Springz>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</Springz>";
  return os;
}
