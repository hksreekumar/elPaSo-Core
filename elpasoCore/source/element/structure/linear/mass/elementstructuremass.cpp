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

#include "elementstructuremass.h"

cElementStructureMass::cElementStructureMass()
    : cElementStructureLinear(1, 3, 0) {
  // empty
}

cElementStructureMass::cElementStructureMass(const cElementStructureMass &other)
    : cElementStructureLinear(other) {
  // empty
}

cElementStructureMass::~cElementStructureMass() {
  // empty
}

std::vector<eKnownDofs> cElementStructureMass::getDofs(void) const {
  std::vector<eKnownDofs> res(3);

  res[0] = disp_x1;
  res[1] = disp_x2;
  res[2] = disp_x3;

  return res;
}

void cElementStructureMass::assembleStiffnessMatrix(cElementMatrix &KM,
                                                    Vec *x = NULL,
                                                    Vec *dx = NULL) {
  // nothing to do here
}

void cElementStructureMass::assembleMassMatrix(cElementMatrix &MM) {
  const PetscReal M = m_Material->getMass();

  MM(0, 0) = M;
  MM(1, 1) = M;
  MM(2, 2) = M;
}

void cElementStructureMass::assembleLoadVector(cElementVector &LV,
                                               cElementMatrix &KM,
                                               Vec *x = NULL, Vec *dx = NULL) {
  // nothing to do here
}

std::ostream &cElementStructureMass::write(std::ostream &os) const {
  os << "mass element" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k = 0; k < getNumberOfNodes(); k++) os << *(m_Nodes[k]);

  return os;
}

std::ostream &cElementStructureMass::writeXml(std::ostream &os) const {
  os << "<Pointmass>";
  os << "<Id>" << getId() << "</Id>";
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</Pointmass>";
  return os;
}
