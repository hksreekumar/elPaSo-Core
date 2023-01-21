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

#include "nodalvaluesvelocity.h"


cNodalValuesVelocity::cNodalValuesVelocity()
{
  m_Velocity[0] = 0.;
  m_Velocity[1] = 0.;
  m_Velocity[2] = 0.;
}


cNodalValuesVelocity::cNodalValuesVelocity(const cNodalValuesVelocity &other) :
  cNodalValues(other)
{
  m_Velocity[0] = other[0];
  m_Velocity[1] = other[1];
  m_Velocity[2] = other[2];

}


cNodalValuesVelocity::~cNodalValuesVelocity()
{
  // empty
}

/*void cNodalValuesVelocity::setPressureValue(PetscScalar Value)
{
  m_P=Value;
} */


// Overload of operator [] needed here and therefore kept 
PetscScalar& cNodalValuesVelocity::operator[](short index)
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_Velocity[index];
}


PetscScalar cNodalValuesVelocity::operator[](short index) const
{
#ifdef PETSC_USE_DEBUG
  if (index<0 || index>2)
    throw cException("bad index", __FILE__, __LINE__);
#endif

  return m_Velocity[index];
}


std::istream& cNodalValuesVelocity::read(std::istream& is)
{
  cId::read(is);
  is >> m_Velocity[0] >> m_Velocity[1] >> m_Velocity[2];

  return is;
}


std::ostream& cNodalValuesVelocity::write(std::ostream &os) const
{
  os << "nodal values velocity" << std::endl;
  os << "  Id : " << getId() << std::endl;
  os << "  Vx : " << m_Velocity[0] << std::endl;
  os << "  Vy : " << m_Velocity[1] << std::endl;
  os << "  Vz : " << m_Velocity[2] << std::endl;

  return os;
}

// output not (yet) needed 
// in femparser not yet test on cNocNodalValues with type "velocity"
std::ostream& cNodalValuesVelocity::writeXml(std::ostream &os) const
{
  os << "<NodalValues Type=\"velocity\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<Vx>" << m_Velocity[0] << "</Vx>";
  os << "<Vy>" << m_Velocity[1] << "</Vy>";
  os << "<Vz>" << m_Velocity[2] << "</Vz>";
  os << "</NodalValues>";

  return os;
}
