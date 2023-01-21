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

#include "elementff.h"

cElementFF::cElementFF(int NumberOfNodes) : cElement(NumberOfNodes) {
  // m_Symmetric = true;
}

cElementFF::cElementFF(const cElementFF &other) : cElement(other) {}

cElementFF::~cElementFF() {}

std::ostream &cElementFF::write(std::ostream &os) const {
  os << "FF (4 or 9 nodes)" << std::endl;
  os << "id = " << getId() << std::endl;
  for (int k = 0; k < getNumberOfNodes(); k++) os << *(m_Nodes[k]) << std::endl;

  return os;
}
