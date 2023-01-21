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

#include "elementstructurespringbcy.h"

cElementStructureSpringBCy::cElementStructureSpringBCy() :
  cElementStructureLinear(1, 1, 0)
{
  // empty
}


cElementStructureSpringBCy::cElementStructureSpringBCy(const cElementStructureSpringBCy &other) :
  cElementStructureLinear(other)
{
  // empty
}


cElementStructureSpringBCy::~cElementStructureSpringBCy()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureSpringBCy::getDofs(void) const
{
  std::vector<eKnownDofs> res(1);

  res[0] = disp_x2;
  return res;
}


void cElementStructureSpringBCy::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  const PetscScalar Cy = m_Material->getCy();
  
  KM(0,0) =  Cy;
}


void cElementStructureSpringBCy::assembleMassMatrix(cElementMatrix &MM)
{
  // nothing to do here
}


void cElementStructureSpringBCy::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  // nothing to do here
}


std::ostream& cElementStructureSpringBCy::write(std::ostream &os) const
{
  os << "spring element bc y" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k=0; k<getNumberOfNodes(); k++)
    os << *(m_Nodes[k]);

  return os;
}


std::ostream& cElementStructureSpringBCy::writeXml(std::ostream &os) const
{
  os << "<SpringBCy>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</SpringBCy>";
  return os;
}
