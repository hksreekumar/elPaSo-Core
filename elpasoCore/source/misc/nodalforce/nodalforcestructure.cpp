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

#include "nodalforcestructure.h"


cNodalForceStructure::cNodalForceStructure()
{
  m_F[0] = 0.;
  m_F[1] = 0.;
  m_F[2] = 0.;
}


cNodalForceStructure::cNodalForceStructure(const cNodalForceStructure &other) :
  cNodalForce(other)
{
  m_F[0] = other[0];
  m_F[1] = other[1];
  m_F[2] = other[2];
}


cNodalForceStructure::~cNodalForceStructure()
{
  // empty
}


PetscReal& cNodalForceStructure::operator[](short index)
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_F[index];
}


PetscReal cNodalForceStructure::operator[](short index) const
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_F[index];
}


std::istream& cNodalForceStructure::read(std::istream& is)
{
  cId::read(is);
  is >> m_F[0] >> m_F[1] >> m_F[2];
  return is;
}


std::ostream& cNodalForceStructure::write(std::ostream &os) const
{
  os << "nodal force structure" << std::endl;
  os << "  Id : " << getId() << std::endl;
  os << "  Fx : " << m_F[0] << std::endl;
  os << "  Fy : " << m_F[1] << std::endl;
  os << "  Fz : " << m_F[2] << std::endl;

  return os;
}


std::ostream& cNodalForceStructure::writeXml(std::ostream &os) const
{
  os << "<NodalForce Type=\"structure\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<Fx>" << m_F[0] << "</Fx>";
  os << "<Fy>" << m_F[1] << "</Fy>";
  os << "<Fz>" << m_F[2] << "</Fz>";
  os << "</NodalForce>";

  return os;
}
