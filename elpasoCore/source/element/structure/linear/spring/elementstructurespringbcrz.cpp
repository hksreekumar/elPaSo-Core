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

#include "elementstructurespringbcrz.h"

cElementStructureSpringBCrz::cElementStructureSpringBCrz() :
  cElementStructureLinear(1, 1, 0)
{
  // empty
}


cElementStructureSpringBCrz::cElementStructureSpringBCrz(const cElementStructureSpringBCrz &other) :
  cElementStructureLinear(other)
{
  // empty
}


cElementStructureSpringBCrz::~cElementStructureSpringBCrz()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureSpringBCrz::getDofs(void) const
{
  std::vector<eKnownDofs> res(1);

  res[0] = disp_w3;
  return res;
}


void cElementStructureSpringBCrz::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  const PetscScalar Crz = m_Material->getCrz();
  
  KM(0,0) =  Crz;
}


void cElementStructureSpringBCrz::assembleMassMatrix(cElementMatrix &MM)
{
  // nothing to do here
}


void cElementStructureSpringBCrz::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  // nothing to do here
}


std::ostream& cElementStructureSpringBCrz::write(std::ostream &os) const
{
  os << "spring element bc rz" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k=0; k<getNumberOfNodes(); k++)
    os << *(m_Nodes[k]);

  return os;
}


std::ostream& cElementStructureSpringBCrz::writeXml(std::ostream &os) const
{
  os << "<SpringBCrz>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</SpringBCrz>";
  return os;
}
