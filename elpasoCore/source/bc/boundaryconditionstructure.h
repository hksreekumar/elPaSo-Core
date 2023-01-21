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

#ifndef INFAM_BOUNDARY_CONDITION_STRUCTURE_H
#define INFAM_BOUNDARY_CONDITION_STRUCTURE_H

#include "boundarycondition.h"

//! @brief nodal boundary condition for nodes of structure elements
//! @author Dirk Clasen
//! @date 24.01.2003
//!
//! The array m_Fixed stores the degrees of freedom in following order
//! @verbatim
//!   m_Fixed[ 0] = u_1
//!   m_Fixed[ 1] = u_2
//!   m_Fixed[ 2] = u_3
//!   m_Fixed[ 3] = w_1
//!   m_Fixed[ 4] = w_2
//!   m_Fixed[ 5] = w_3
//!   m_Fixed[ 6] = p_0 (porepressure)
//!   m_Fixed[ 7] = p_1 (porepressure)
//!   m_Fixed[ 8] = p_2 (porepressure)
//!   m_Fixed[ 9] = p_3 (porepressure)
//!   m_Fixed[10] = xd3,
//!   m_Fixed[11] = wd1,
//!   m_Fixed[12] = wd2,
//!   m_Fixed[13] = dwdx, Kirchhoff
//!   m_Fixed[14] = dwdy, Kirchhoff
//!   m_Fixed[15] = dwdxy Kirchhoff
//!   m_Fixed[16] = fluid
//!   m_Fixed[17] = disp_z_1
//!   m_Fixed[18] = disp_z_3
//!   m_Fixed[19] = disp_x1_2
//!   m_Fixed[20] = disp_x2_2
//! @endverbatim
//! \f[ p = \sum_{i=0}^n p_i x^i \f]
//! The indices of the array m_Value are sorted in the same manner.
class cBoundaryConditionStructure
    : public cBoundaryCondition,
      private tCounter<cBoundaryConditionStructure> {
 private:
  enum {
    eNumberOfDofs =
        21  ///< max. number of degrees of freedom at structureelements
  };

  std::bitset<eNumberOfDofs>
      m_Fixed;  ///< tells us, if a degree of freedom is fixed
  std::vector<PetscReal>
      m_Value;  ///< prescribed value for the associated degree of freedom

 public:
  //! @brief Constructor
  cBoundaryConditionStructure();
  //! @brief Copy constructor
  cBoundaryConditionStructure(const cBoundaryConditionStructure &other);
  //! @brief Destructor
  ~cBoundaryConditionStructure();

  //! @brief returns the number of instances of the current object
  static size_t howMany(void) {
    return tCounter<cBoundaryConditionStructure>::howMany();
  }

  //! @brief check if a degree of freedom is fixed
  //! @param dof denomination of the degree of freedom
  //! @see testboundaryconditionstructure.h
  bool checkIfFixed(eKnownDofs dof) const;
  //! @brief check if a degree of freedom is fixed
  //! @param dof id of the degree of freedom
  //! @see testboundaryconditionstructure.h
  bool checkIfFixed(int dof) const { return m_Fixed.test(dof); };

  //! @brief returns the prescriebed value of a specific degree of freedom
  //! @param dof denomination of the degree of freedom
  //! @see testboundaryconditionstructure.h
  PetscScalar getPrescribedValue(eKnownDofs dof, cPoint *curPt,
                                 PetscReal omega) const;

  //! @brief tell a degree of freedom that he is fixed
  //! @param dof index of the degree of freedom to fix
  //! @see testboundaryconditionstructure.h
  inline void fixDof(int dof) { m_Fixed.set(dof, 1); }

  //! @brief prescribe a value for a specific degree of freedom
  //! @param dof index of the degree of freedom
  //! @param value value to be prescribed
  //! @see testboundaryconditionstructure.h
  inline void setPrescribedValue(int dof, const PetscReal &value) {
    m_Value[dof] = value;
  }

  //! @brief returns the whole array - used by the copy constructor
  //! @see testboundaryconditionstructure.h
  inline std::bitset<eNumberOfDofs> getFixedArray(void) const {
    return m_Fixed;
  }

  //! @brief returns the whole array - used by the copy constructor
  //! @see testboundaryconditionstructure.h
  inline std::vector<PetscReal> getValueArray(void) const { return m_Value; }

  /*BEGIN_NO_COVERAGE*/
  //! @brief writes this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief reads data to this object
  //! @param is inputstream
  //! @return modified inputstream
  std::istream &read(std::istream &is);

  //! @brief writes this object in XML
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
  /*END_NO_COVERAGE*/
};

//! @brief overloaded output operator
inline std::ostream &operator<<(std::ostream &os,
                                const cBoundaryConditionStructure &other) {
  return other.write(os);
}

//! @brief overloaded input operator
inline std::istream &operator>>(std::istream &is,
                                cBoundaryConditionStructure &other) {
  return other.read(is);
}

#endif
