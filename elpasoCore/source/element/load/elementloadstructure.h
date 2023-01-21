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

#ifndef INFAM_ELEMENT_LOAD_STRUCTURE_H
#define INFAM_ELEMENT_LOAD_STRUCTURE_H

#include "elementload.h"

//! @brief base class for element loads of structure elements
//! @author Dirk Clasen
//! @date 31.05.2005
class cElementLoadStructure : public cElementLoad,
                              private tCounter<cElementLoadStructure> {
 public:
  //! @brief Constructor
  cElementLoadStructure();
  //! @brief Copy constructor
  cElementLoadStructure(const cElementLoadStructure &other);
  //! @brief Destructor
  virtual ~cElementLoadStructure();

  //! @brief return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementLoadStructure>::howMany();
  }

  //! @brief return a single component of the forcevector
  //! @param NodeId id of node to whom load will be applied
  //! @param i direction in global coordinate system of the sought-for component
  //! of the force
  //! @return entry of force vector
  virtual PetscScalar getForceComponent(PetscInt NodeId, int i) const = 0;

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief read data from a stream to this object
  //! @param is inputstream
  //! @return modified inputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! @brief write this object to a XML stream
  //! @param os der Ausgabestream
  //! @return der modifizierte Ausgabestream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementLoadStructure &other) {
  return other.write(os);
}

//! @brief overloaded inputoperator
inline std::istream &operator>>(std::istream &is,
                                cElementLoadStructure &other) {
  return other.read(is);
}

#endif
