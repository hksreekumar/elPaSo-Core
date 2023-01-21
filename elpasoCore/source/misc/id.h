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

#ifndef INFAM_ID_H
#define INFAM_ID_H

#include "mytypes.h"

//! @brief stores the identifier of an object of the mesh
//! @author Dirk Clasen
//! @date 04.08.2005
class cId {
 private:
  PetscInt m_Id;    ///< identifier of the object: index form x to y
  PetscInt m_Id0n;  ///< identifier of the object: index from 0 to n

 public:
  cId();
  cId(const cId &other);
  virtual ~cId();

  //! assigns a new id to an object
  //! @param NewId the new id
  void setId(PetscInt NewId);

  //! assigns a new id to an object
  //! @param NewId0n the new id0n
  void setId0n(PetscInt NewId0n);

  //! returns the id of an object
  inline PetscInt getId(void) const { return m_Id; }

  //! returns the id0n of an object
  inline PetscInt getId0n(void) const { return m_Id0n; }

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified intputstream
  std::istream &read(std::istream &is);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cId &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cId &other) {
  return other.write(os);
}

#endif
