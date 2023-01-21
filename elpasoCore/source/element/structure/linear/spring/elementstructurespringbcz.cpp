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

#include "elementstructurespringbcz.h"

cElementStructureSpringBCz::cElementStructureSpringBCz() :
  cElementStructureLinear(1, 1, 0)
{
  // empty
}


cElementStructureSpringBCz::cElementStructureSpringBCz(const cElementStructureSpringBCz &other) :
  cElementStructureLinear(other)
{
  // empty
}


cElementStructureSpringBCz::~cElementStructureSpringBCz()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureSpringBCz::getDofs(void) const
{
  std::vector<eKnownDofs> res(1);

  res[0] = disp_x3;
  return res;
}


void cElementStructureSpringBCz::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  const PetscScalar Cz = m_Material->getCz();
  
  KM(0,0) =  Cz;
}


void cElementStructureSpringBCz::assembleMassMatrix(cElementMatrix &MM)
{
  // nothing to do here
}


void cElementStructureSpringBCz::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  // nothing to do here
}


std::ostream& cElementStructureSpringBCz::write(std::ostream &os) const
{
  os << "spring element bc z" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k=0; k<getNumberOfNodes(); k++)
    os << *(m_Nodes[k]);

  return os;
}


std::ostream& cElementStructureSpringBCz::writeXml(std::ostream &os) const
{
  os << "<SpringBCz>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</SpringBCz>";
  return os;
}
