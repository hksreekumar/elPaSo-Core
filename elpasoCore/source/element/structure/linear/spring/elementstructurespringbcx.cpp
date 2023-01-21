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

#include "elementstructurespringbcx.h"

cElementStructureSpringBCx::cElementStructureSpringBCx() :
  cElementStructureLinear(1, 1, 0)
{
  // empty
}


cElementStructureSpringBCx::cElementStructureSpringBCx(const cElementStructureSpringBCx &other) :
  cElementStructureLinear(other)
{
  // empty
}


cElementStructureSpringBCx::~cElementStructureSpringBCx()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureSpringBCx::getDofs(void) const
{
  std::vector<eKnownDofs> res(1);

  res[0] = disp_x1;
  return res;
}


void cElementStructureSpringBCx::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  const PetscScalar Cx = m_Material->getCx();
  
  KM(0,0) =  Cx;
}


void cElementStructureSpringBCx::assembleMassMatrix(cElementMatrix &MM)
{
  // nothing to do here
}


void cElementStructureSpringBCx::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  // nothing to do here
}


std::ostream& cElementStructureSpringBCx::write(std::ostream &os) const
{
  os << "spring element bc x" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k=0; k<getNumberOfNodes(); k++)
    os << *(m_Nodes[k]);

  return os;
}


std::ostream& cElementStructureSpringBCx::writeXml(std::ostream &os) const
{
  os << "<SpringBCx>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</SpringBCx>";
  return os;
}
