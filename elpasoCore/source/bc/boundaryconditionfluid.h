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

#ifndef INFAM_BOUNDARY_CONDITION_FLUID_H
#define INFAM_BOUNDARY_CONDITION_FLUID_H

#include "boundarycondition.h"

//! @brief boundary condition for fluid nodes
//! @author Dirk Clasen
//! @date 03.02.2005
class cBoundaryConditionFluid : public cBoundaryCondition,
                                private tCounter<cBoundaryConditionFluid> {
 private:
  PetscReal m_Pressure;  ///< prescribed fluid pressure

 public:
  //! @brief Constructor
  cBoundaryConditionFluid();
  //! @brief Copy constructor
  cBoundaryConditionFluid(const cBoundaryConditionFluid &other);
  //! @brief Destructor
  ~cBoundaryConditionFluid();

  //! @brief read the instance counter of this object type
  static size_t howMany(void) {
    return tCounter<cBoundaryConditionFluid>::howMany();
  }

  //! @brief checks if a degree of freedom is fixed
  //! @param dof enum dof identifier
  //! @return true if the dof is fixed; else false
  //! @see testboundaryconditionfluid.h
  bool checkIfFixed(eKnownDofs dof) const;
  //! @brief checks if a degree of freedom is fixed
  //! @param dof id
  //! @return true if the dof is fixed; else false
  //! @see testboundaryconditionfluid.h
  bool checkIfFixed(int dof) const;

  //! @brief reads the prescribed value for a fixed degree of freedom
  //! @see testboundaryconditionfluid.h
  virtual PetscScalar getPrescribedValue(eKnownDofs dof, cPoint *curPt,
                                         PetscReal omega) const;

  //! @brief set the value for the prescribed pressure
  //! @see testboundaryconditionfluid.h
  void setPressure(const PetscReal &pressure) { m_Pressure = pressure; }

  //! @brief returns the prescribed fluid pressure
  //! @see testboundaryconditionfluid.h
  PetscReal getPressure(void) const { return m_Pressure; }

  //! @brief output of this to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const;

  //! @brief reads this object from a stream
  //! @param is inputstream
  //! @return modified inputstream
  virtual std::istream &read(std::istream &is);

  //! @brief writes this object in XML format
  //! This function will be used to export the mesh to an
  //! inputfile after looking for coupling elements
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const;
};

//! @brief overloaded output operator
inline std::ostream &operator<<(std::ostream &os,
                                const cBoundaryConditionFluid &other) {
  return other.write(os);
}

//! @brief overloaded inputoperator
inline std::istream &operator>>(std::istream &is,
                                cBoundaryConditionFluid &other) {
  return other.read(is);
}

#endif
