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

#include "id.h"

cId::cId(void)
{
  setId(0);
  setId0n(0);
}

cId::cId(const cId &other)
{
  setId(other.getId());
  setId0n(other.getId0n());
}

cId::~cId()
{
  // leer
}


void cId::setId(PetscInt NewId)
{
#ifdef PETSC_USE_DEBUG
  if (NewId < 0)
    throw cException("Id kleiner Null", __FILE__, __LINE__);
#endif

  m_Id = NewId;
}

void cId::setId0n(PetscInt NewId0n)
{
#ifdef PETSC_USE_DEBUG
  if (NewId0n < 0)
    throw cException("Id0n kleiner Null", __FILE__, __LINE__);
#endif

  m_Id0n = NewId0n;
}

std::istream& cId::read(std::istream &is)
{
  is >> m_Id;
  return is;
}

std::ostream& cId::write(std::ostream &os) const
{
  os << std::setw(5) << m_Id;
  return os;
}

