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

#include "nodalpressurestructure.h"


cNodalPressureStructure::cNodalPressureStructure()
{
  m_P= (0.,0.);
}


cNodalPressureStructure::cNodalPressureStructure(const cNodalPressureStructure &other) :
  cNodalPressure(other)
{
  m_P = other.m_P;
}


cNodalPressureStructure::~cNodalPressureStructure()
{
  // empty
}

void cNodalPressureStructure::setPressureValue(PetscScalar Value)
{
  m_P=Value;
}


// Overload of operator [] not really needed here but kept for possible compatibility issues
PetscScalar& cNodalPressureStructure::operator[](short index)
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_P;
}


PetscScalar cNodalPressureStructure::operator[](short index) const
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_P;
}


std::istream& cNodalPressureStructure::read(std::istream& is)
{
  cId::read(is);
  is >> m_P;

  return is;
}


std::ostream& cNodalPressureStructure::write(std::ostream &os) const
{
  os << "nodal pressure structure" << std::endl;
  os << "  Id : " << getId() << std::endl;
  os << "  P  : " << m_P << std::endl;

  return os;
}


std::ostream& cNodalPressureStructure::writeXml(std::ostream &os) const
{
  os << "<NodalPressure Type=\"structure\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<P>" << m_P << "</P>";
  os << "</NodalPressure>";

  return os;
}
