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

#ifndef INFAM_NODALVALUES_H
#define INFAM_NODALVALUES_H

#include "../../basics/exceptions/exceptions.h"
#include "../counter.h"
#include "../id.h"

//! @brief class for nodal values
//! @author Silja Beck
//! @date 09.09.2011
//! Class cNodalValues is build the same way as cNodalForces
class cNodalValues : public cId,
                     public virtual cLogging,
                     private tCounter<cNodalValues> {
 private:
 public:
  cNodalValues();
  cNodalValues(const cNodalValues &other);
  virtual ~cNodalValues();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cNodalValues>::howMany(); }

  // PetscScalar getPressureComponent();

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cNodalValues &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cNodalValues &other) {
  return other.write(os);
}

#endif
