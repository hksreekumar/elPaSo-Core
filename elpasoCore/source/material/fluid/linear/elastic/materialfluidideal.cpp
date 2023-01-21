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

#include "materialfluidideal.h"

/*BEGIN_NO_COVERAGE*/
cMaterialFluidIdeal::cMaterialFluidIdeal()
{
  // empty
}


cMaterialFluidIdeal::cMaterialFluidIdeal(const cMaterialFluidIdeal &other) :
  cMaterialFluidLin(other)
{
  // empty
}

cMaterialFluidIdeal::~cMaterialFluidIdeal()
{
  // empty
}


std::istream& cMaterialFluidIdeal::read(std::istream &is)
{
  cId::read(is);
  is >> m_Cf >> m_Rho >> m_T;

  return is;
}


std::ostream& cMaterialFluidIdeal::write(std::ostream &os) const
{
  os << "Fluid-Material (" << getIdentifier() << ")" << std::endl;
  os << "  Id  = " << getId() << std::endl;
  os << "  c   = " << getCf() << std::endl;
  os << "  rho = " << getRho() << std::endl;
  os << "  t   = " << getT() << std::endl;
  return os;
}


std::ostream& cMaterialFluidIdeal::writeXml(std::ostream &os) const
{
  os << "<Material Type=\"fluid\" Name=\"" << getIdentifier() << "\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<cf>" << getCf() << "</cf>";
  os << "<rhof>" << getRho() << "</rhof>";
  os << "<t>" << getT() << "</t>";
  os << "</Material>";

  return os;
}
/*END_NO_COVERAGE*/