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

#ifndef INFAM_NODALVALUES_VELOCITY_H
#define INFAM_NODALVALUES_VELOCITY_H

#include "nodalvalues.h"

//! @brief class for nodal values velocity
class cNodalValuesVelocity : public cNodalValues,
                             private tCounter<cNodalValuesVelocity> {
 private:
  PetscScalar m_Velocity[3];  ///< velocity values (3d) for one node from fluid
                              ///< flow (NodalValuesFF)

 public:
  cNodalValuesVelocity();
  cNodalValuesVelocity(const cNodalValuesVelocity &other);
  ~cNodalValuesVelocity();

  //! return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cNodalValuesVelocity>::howMany();
  }

  //! function to get velocity value at the node in dimension d
  inline PetscScalar getVelocityComponent(int k) const { return m_Velocity[k]; }

  //! function to set velocity value at the node
  // void setVelocityValue(PetscScalar Value);

  // Overload of operator [] not really needed here but kept for possible
  // compatibility issues
  //! access one component of the force (read and write)
  PetscScalar &operator[](short index);

  //! access one component of the force (read-only)
  PetscScalar operator[](short index) const;

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cNodalValuesVelocity &other) {
  return other.write(os);
}

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cNodalValuesVelocity &other) {
  return other.read(is);
}

#endif
