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

#ifndef INFAM_NODALFORCE_FLUID_H
#define INFAM_NODALFORCE_FLUID_H

#include "nodalforce.h"

//! @brief Base class for nodal force fluid domain
class cNodalForceFluid : public cNodalForce,
                         private tCounter<cNodalForceFluid> {
 private:
  PetscReal m_Source;  ///< point sound source

 public:
  cNodalForceFluid();
  cNodalForceFluid(const cNodalForceFluid &other);
  ~cNodalForceFluid();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cNodalForceFluid>::howMany(); }

  inline PetscReal getSourceValue(void) const { return m_Source; }

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
                                const cNodalForceFluid &other) {
  return other.write(os);
}

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cNodalForceFluid &other) {
  return other.read(is);
}

#endif
