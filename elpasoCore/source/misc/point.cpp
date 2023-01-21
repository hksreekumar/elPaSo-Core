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

#include "point.h"

cPoint::cPoint(void)
{
  m_X[0] = 0.;
  m_X[1] = 0.;
  m_X[2] = 0.;
}

cPoint::cPoint(const cPoint &other)
{
  m_X[0] = other[0];
  m_X[1] = other[1];
  m_X[2] = other[2];
}

/*BEGIN_NO_COVERAGE*/
cPoint::~cPoint()
{
  // empty
}
/*END_NO_COVERAGE*/

PetscReal cPoint::operator[](int i) const
{
#ifdef PETSC_USE_DEBUG
  if (i<0 || i>2)
    throw cException("invalid index", __FILE__, __LINE__);
#endif

  return m_X[i];
}


PetscReal& cPoint::operator[](int i)
{
#ifdef PETSC_USE_DEBUG
  if (i<0 || i>2)
    throw cException("invalid index", __FILE__, __LINE__);
#endif

  return m_X[i];
}

PetscReal cPoint::getComponent(int i)
{
#ifdef PETSC_USE_DEBUG
  if (i<0 || i>2)
    throw cException("invalid index", __FILE__, __LINE__);
#endif

  return m_X[i];
}

PetscReal cPoint::distance(const cPoint &other) const
{
  PetscReal r = 0.;

  for (int k=0; k<3; k++)
    r += (m_X[k] - other[k]) * (m_X[k] - other[k]);

  return std::sqrt(r);
}


PetscReal cPoint::distance2(const cPoint &other) const
{
  PetscReal r = 0.;

  for (int k=0; k<3; k++)
    r += (m_X[k] - other[k]) * (m_X[k] - other[k]);

  return r;
}

/*BEGIN_NO_COVERAGE*/
std::istream& cPoint::read(std::istream &is)
{
  is >> m_X[0] >> m_X[1] >> m_X[2];
  return is;
}

std::ostream& cPoint::write(std::ostream &os) const
{
  os.setf(std::ios::scientific);
  os << std::setw(15) << m_X[0];
  os << std::setw(15) << m_X[1];
  os << std::setw(15) << m_X[2];
  os.unsetf(std::ios::scientific);

  return os;
}
/*END_NO_COVERAGE*/
