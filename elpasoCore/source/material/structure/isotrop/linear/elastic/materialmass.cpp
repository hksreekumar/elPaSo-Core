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

#include "materialmass.h"

cMaterialMass::cMaterialMass(void) :
  cMaterialStructureIsoLin()
{
  m_M = 0.;
}


cMaterialMass::cMaterialMass(const cMaterialMass &other) :
  cMaterialStructureIsoLin(other)
{
  m_M = other.getMass();
}


cMaterialMass::~ cMaterialMass()
{
  // leer
}


std::istream& cMaterialMass::read(std::istream &is)
{
  cId::read(is);
  is >> m_M;

  return is;
}


std::ostream& cMaterialMass::write(std::ostream &os) const
{
  os << "Mass material (" << getIdentifier() << ")" << std::endl;
  os << "Id..: " << getId() << std::endl;
  os << "M..: " << getMass() << std::endl;

  return os;
}


std::ostream& cMaterialMass::writeXml(std::ostream &os) const
{
  os << "<Material Type=\"pointmass\" Name=\"" << getIdentifier() << "\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<M>" << getCx() << "</M>";
  os << "</Material>";

  return os;
}
