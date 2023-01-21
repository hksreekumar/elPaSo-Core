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

#ifndef INFAM_ELEMENT_LOAD_STRUCTURE_CONST_H
#define INFAM_ELEMENT_LOAD_STRUCTURE_CONST_H

#include "elementloadstructure.h"

//! @brief base class for distributed element loads
//! @author Dirk Clasen
//! @date 05.10.2005
class cElementLoadStructureConst
    : public cElementLoadStructure,
      private tCounter<cElementLoadStructureConst> {
 private:
  PetscScalar m_F[6];
  PetscReal m_Phase;

 public:
  //! @brief Constructor
  cElementLoadStructureConst();
  //! @brief Copy constructor
  cElementLoadStructureConst(const cElementLoadStructureConst &other);
  //! @brief Destructor
  ~cElementLoadStructureConst();

  //! @brief return the number of instances of this class
  static size_t howMany(void) {
    return tCounter<cElementLoadStructureConst>::howMany();
  }

  //! @brief return a single component of the forcevector
  //! @param NodeId id of node to whom load will be applied
  //! @param i direction in global coordinate system of the sought-for component
  //! of the force
  //! @return entry of force vector
  PetscScalar getForceComponent(PetscInt NodeId, int i) const { return m_F[i]; }

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief read data from a stream to this object
  //! @param is inputstream
  //! @return modified inputstream
  std::istream &read(std::istream &is);

  //! @brief write this object to a XML stream
  //! @param os der Ausgabestream
  //! @return der modifizierte Ausgabestream
  std::ostream &writeXml(std::ostream &os) const;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementLoadStructureConst &other) {
  return other.write(os);
}

//! @brief overloaded inputoperator
inline std::istream &operator>>(std::istream &is,
                                cElementLoadStructureConst &other) {
  return other.read(is);
}

#endif
