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

#include "elementstructurespringbc.h"

cElementStructureSpringBC::cElementStructureSpringBC() :
  cElementStructureLinear(1, 6, 0)
{
  // empty
}


cElementStructureSpringBC::cElementStructureSpringBC(const cElementStructureSpringBC &other) :
  cElementStructureLinear(other)
{
  // empty
}


cElementStructureSpringBC::~cElementStructureSpringBC()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureSpringBC::getDofs(void) const
{
  std::vector<eKnownDofs> res(6);

  res[0] = disp_x1;
  res[1] = disp_x2;
  res[2] = disp_x3;
  res[3] = disp_w1;
  res[4] = disp_w2;
  res[5] = disp_w3;
  return res;
}


void cElementStructureSpringBC::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  const PetscScalar Cx = m_Material->getCx();
  const PetscScalar Cy = m_Material->getCy();
  const PetscScalar Cz = m_Material->getCz();
  const PetscScalar Crx = m_Material->getCrx();
  const PetscScalar Cry = m_Material->getCry();
  const PetscScalar Crz = m_Material->getCrz();
  KM(0,0) =  Cx;
          KM(1,1) =  Cy;
                  KM(2,2) =  Cz;
                          KM(3,3) =  Crx;
                                  KM(4,4) =  Cry;
                                          KM(5,5) =  Crz;
}


void cElementStructureSpringBC::assembleMassMatrix(cElementMatrix &MM)
{
  // nothing to do here
}


void cElementStructureSpringBC::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  // nothing to do here
}


std::ostream& cElementStructureSpringBC::write(std::ostream &os) const
{
  os << "spring element bc" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k=0; k<getNumberOfNodes(); k++)
    os << *(m_Nodes[k]);

  return os;
}


std::ostream& cElementStructureSpringBC::writeXml(std::ostream &os) const
{
  os << "<SpringBC>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</SpringBC>";
  return os;
}
