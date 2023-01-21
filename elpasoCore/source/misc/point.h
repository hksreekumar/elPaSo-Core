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

#ifndef INFAM_POINT_H
#define INFAM_POINT_H

#include <cmath>

#include "mytypes.h"

//! @brief 3d point within cartesian coordinates
//! @author Dirk Clasen
//! @date 26.07.2005
class cPoint {
 private:
  PetscReal m_X[3];  ///< components x_i of the point

 public:
  //! @see testpoint.h
  cPoint(void);
  //! @see testpoint.h
  cPoint(const cPoint &other);
  virtual ~cPoint();

  //! access one component (read-only)
  //! @see testpoint.h
  PetscReal operator[](int i) const;

  //! access one component of the point (read/write)
  //! @see testpoint.h
  PetscReal &operator[](int i);

  //! return component i of the point
  //! @see testpoint.h
  PetscReal getComponent(int i);

  //! computes the distance between this and other
  //! @see testpoint.h
  PetscReal distance(const cPoint &other) const;

  //! computes the squared distance between this and other
  //! @see testpoint.h
  PetscReal distance2(const cPoint &other) const;

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cPoint &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cPoint &other) {
  return other.write(os);
}

#endif
