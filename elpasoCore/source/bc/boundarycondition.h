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

#ifndef INFAM_BOUNDARY_CONDITION_H
#define INFAM_BOUNDARY_CONDITION_H

#include "../fedata/dof.h"
#include "../misc/counter.h"
#include "../misc/id.h"
#include "../misc/point.h"

//! @brief virtual baseclass of all nodal boundaryconditions
//! @author Dirk Clasen
//! @date 02.06.2005
class cBoundaryCondition : public cId, private tCounter<cBoundaryCondition> {
 protected:
  std::string m_Identifier;  ///< identifier of the boundary condition (as
                             ///< specified in Patran)

 public:
  //! @brief Constructor
  cBoundaryCondition();
  //! @brief Copy constructor
  cBoundaryCondition(const cBoundaryCondition &other);
  //! @brief Destructor
  virtual ~cBoundaryCondition();

  //! @brief returns the number of instances of this objecttype
  static size_t howMany(void) {
    return tCounter<cBoundaryCondition>::howMany();
  }

  //! @brief checks if a degree of freedom is fixed
  virtual bool checkIfFixed(eKnownDofs dof) const = 0;
  virtual bool checkIfFixed(int dof) const = 0;

  //! @brief reads the prescribed value for a fixed degree of freedom
  virtual PetscScalar getPrescribedValue(eKnownDofs dof, cPoint *curPt,
                                         PetscReal omega) const = 0;

  //! @brief return the identifier of this boundary condition
  //! @see testboundaryconditionfluid.h, testboundaryconditionstructure.h
  inline std::string getIdentifier(void) const { return m_Identifier; }

  //! @brief specify a new identifier for this boundary condition
  inline void setIdentifier(const std::string &name) { m_Identifier = name; }

  //! @brief writes this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief reads a single object of a stream
  //! @param is inputstream
  //! @return modified inputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! @brief writes this object in XML format
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cBoundaryCondition &other) {
  return other.write(os);
}

//! overloaded outputoperator
inline std::istream &operator>>(std::istream &is, cBoundaryCondition &other) {
  return other.read(is);
}

#endif
