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

#include "elementstructurespring.h"

cElementStructureSpring::cElementStructureSpring() :
  cElementStructureLinear(2, 6, 0)
{
  // empty
}


cElementStructureSpring::cElementStructureSpring(const cElementStructureSpring &other) :
  cElementStructureLinear(other)
{
  // empty
}


cElementStructureSpring::~cElementStructureSpring()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureSpring::getDofs(void) const
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


void cElementStructureSpring::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  const PetscScalar Cx = m_Material->getCx();
  const PetscScalar Cy = m_Material->getCy();
  const PetscScalar Cz = m_Material->getCz();
  const PetscScalar Crx = m_Material->getCrx();
  const PetscScalar Cry = m_Material->getCry();
  const PetscScalar Crz = m_Material->getCrz();
  KM(0,0) =  Cx;                              KM(0,6) =  -Cx; 	
          KM(1,1) =  Cy;                              KM(1,7) =  -Cy;
                  KM(2,2) =  Cz;                              KM(2,8) =  -Cz; 
                          KM(3,3) =  Crx;                              KM(3,9) =  -Crx; 
                                  KM(4,4) =  Cry;                              KM(4,10) =  -Cry; 
                                          KM(5,5) =  Crz;                                KM(5,11) =  -Crz; 
  KM(6,0) =  -Cx;                             KM(6,6) =  Cx; 	
          KM(7,1) =  -Cy;                             KM(7,7) =  Cy;
                  KM(8,2) =  -Cz;                             KM(8,8) =  Cz; 
                          KM(9,3) =  -Crx;                             KM(9,9) =  Crx; 
                                  KM(10,4) =  -Cry;                            KM(10,10) =  Cry;
                                            KM(11,5) =  -Crz;                            KM(11,11) =  Crz;     
  //std::cout << KM << std::endl;
}



void cElementStructureSpring::assembleMassMatrix(cElementMatrix &MM)
{
  // nothing to do here
}


void cElementStructureSpring::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
   // nothing to do here
}


std::ostream& cElementStructureSpring::write(std::ostream &os) const
{
  os << "spring element" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k=0; k<getNumberOfNodes(); k++)
    os << *(m_Nodes[k]);

  return os;
}


std::ostream& cElementStructureSpring::writeXml(std::ostream &os) const
{
  os << "<Spring>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</Spring>";
  return os;
}
