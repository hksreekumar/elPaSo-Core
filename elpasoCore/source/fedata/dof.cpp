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

#include "dof.h"

cDegreesOfFreedom::cDegreesOfFreedom()
    : m_GlobalRows(cstNumberOfKnownDofs, -1) {
  m_DofActive.reset();  // set all bits to zero
}

cDegreesOfFreedom::cDegreesOfFreedom(const cDegreesOfFreedom &other)
    : m_DofActive(other.getAllDofFlags()),
      m_GlobalRows(other.getAllGlobalRows()) {
  // empty
}

/*BEGIN_NO_COVERAGE*/
cDegreesOfFreedom::~cDegreesOfFreedom() {
  // empty
}

std::istream &cDegreesOfFreedom::read(std::istream &is) { return is; }

std::ostream &cDegreesOfFreedom::write(std::ostream &os) const {
  for (int k = 0; k < cstNumberOfKnownDofs; k++)
    os << std::setw(7) << "(" << checkIfActive(k) << ":" << getGlobalRow(k)
       << ")";

  return os;
}
/*END_NO_COVERAGE*/