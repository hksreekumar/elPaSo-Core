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

#include "nodalmomentstructure.h"


cNodalMomentStructure::cNodalMomentStructure()
{
  m_M[0] = 0.;
  m_M[1] = 0.;
  m_M[2] = 0.;
}


cNodalMomentStructure::cNodalMomentStructure(const cNodalMomentStructure &other) :
  cNodalMoment(other)
{
  m_M[0] = other[0];
  m_M[1] = other[1];
  m_M[2] = other[2];
}


cNodalMomentStructure::~cNodalMomentStructure()
{
  // empty
}


PetscReal& cNodalMomentStructure::operator[](short index)
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_M[index];
}


PetscReal cNodalMomentStructure::operator[](short index) const
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_M[index];
}


std::istream& cNodalMomentStructure::read(std::istream& is)
{
  cId::read(is);
  is >> m_M[0] >> m_M[1] >> m_M[2];
  return is;
}


std::ostream& cNodalMomentStructure::write(std::ostream &os) const
{
  os << "nodal moment structure" << std::endl;
  os << "  Id : " << getId() << std::endl;
  os << "  Mx : " << m_M[0] << std::endl;
  os << "  My : " << m_M[1] << std::endl;
  os << "  Mz : " << m_M[2] << std::endl;

  return os;
}


std::ostream& cNodalMomentStructure::writeXml(std::ostream &os) const
{
  os << "<NodalMoment Type=\"structure\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<Mx>" << m_M[0] << "</Mx>";
  os << "<My>" << m_M[1] << "</My>";
  os << "<Mz>" << m_M[2] << "</Mz>";
  os << "</NodalMoment>";

  return os;
}
