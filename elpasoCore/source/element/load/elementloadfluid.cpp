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

#include "elementloadfluid.h"

cElementLoadFluid::cElementLoadFluid() { m_Source = 0.; }

cElementLoadFluid::cElementLoadFluid(const cElementLoadFluid &other)
    : cElementLoad(other) {
  m_Source = other.getSourceValue();
}

cElementLoadFluid::~cElementLoadFluid() {
  // empty
}

std::ostream &cElementLoadFluid::write(std::ostream &os) const {
  os << "Fluid-load" << std::endl;
  os << "  Id = " << getId() << std::endl;
  os << "  P  = " << getSourceValue() << std::endl;

  return os;
}

std::istream &cElementLoadFluid::read(std::istream &is) {
  cId::read(is);
  is >> m_Source;
  return is;
}

std::ostream &cElementLoadFluid::writeXml(std::ostream &os) const {
  os << "<ElemLoad Type=\"fluid\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<P>" << getSourceValue() << "</P>";
  os << "</ElemLoad>";

  return os;
}
