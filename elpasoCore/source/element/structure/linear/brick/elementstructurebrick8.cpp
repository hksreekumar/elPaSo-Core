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

#include "elementstructurebrick8.h"

cElementStructureBrick8::cElementStructureBrick8()
    : cElementStructureBrick(8, 2) {
  // empty
}

cElementStructureBrick8::cElementStructureBrick8(
    const cElementStructureBrick8 &other)
    : cElementStructureBrick(other) {
  // empty
}

cElementStructureBrick8::~cElementStructureBrick8() {
  // empty
}

void cElementStructureBrick8::getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
  // -- local coordinates of element's nodes
  xyz(0, 0) = -1.;
  xyz(0, 1) = -1.;
  xyz(0, 2) = -1.;  // node 1
  xyz(1, 0) = 1.;
  xyz(1, 1) = -1.;
  xyz(1, 2) = -1.;  // node 2 ...
  xyz(2, 0) = 1.;
  xyz(2, 1) = 1.;
  xyz(2, 2) = -1.;
  xyz(3, 0) = -1.;
  xyz(3, 1) = 1.;
  xyz(3, 2) = -1.;
  xyz(4, 0) = -1.;
  xyz(4, 1) = -1.;
  xyz(4, 2) = 1.;
  xyz(5, 0) = 1.;
  xyz(5, 1) = -1.;
  xyz(5, 2) = 1.;
  xyz(6, 0) = 1.;
  xyz(6, 1) = 1.;
  xyz(6, 2) = 1.;
  xyz(7, 0) = -1.;
  xyz(7, 1) = 1.;
  xyz(7, 2) = 1.;  // node 8
}

std::ostream &cElementStructureBrick8::write(std::ostream &os) const {
  os << "Brick8 solid element" << std::endl;
  os << "Id = " << getId() << std::endl;
  os << *m_Material << std::endl;

  for (int k = 0; k < getNumberOfNodes(); k++) os << *(m_Nodes[k]) << std::endl;

  return os;
}

std::ostream &cElementStructureBrick8::writeXml(std::ostream &os) const {
  os << "<Brick8>";
  os << "<Id>" << getId() << "</Id>";
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</Brick8>";
  return os;
}
