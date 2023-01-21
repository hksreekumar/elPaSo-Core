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

#ifndef INFAM_DOF_H
#define INFAM_DOF_H

#include "../misc/mytypes.h"

// holds definitions (constants and enums for DOFs)
#include "basics_dof.h"

//! @brief holds the degrees of freedom of a single node
//! This class governs the degrees of freedom at a single nodes.
//! The bitset is used to mark active degrees of freedom by setting
//! the corresponding bit to one. The array m_GlobalRows stores the
//! information about the global row/column of the degree of freedom
//! within the systemmatrices and loadvectors. This information is
//! used for assembling the matrices/vectors.
//! By default, the values of m_GlobalRows are -1. In this case,
//! PETSc ignores them when trying to add them to the matrices/vectors.
//! @author Dirk Clasen
//! @date 22.06.2005
class cDegreesOfFreedom {
 private:
  std::bitset<cstNumberOfKnownDofs>
      m_DofActive;  ///< array of degrees of freedom (active=1, inactive=0)
  std::vector<PetscInt>
      m_GlobalRows;  ///< global rows/columns of a specific degree of freedom

 public:
  cDegreesOfFreedom();
  cDegreesOfFreedom(const cDegreesOfFreedom &other);
  virtual ~cDegreesOfFreedom();

  //! tell a degree of freedom that it exists
  //! @note unit-tested
  inline void activateDof(int dof) { m_DofActive.set(dof, 1); }

  //! tell a degree of freedom that it not exists
  //! @note unit-tested
  inline void deactivateDof(int dof) { m_DofActive.set(dof, 0); }

  //! check if a degree of freedom is marked as active
  //! @param pos degree of freedom to check
  //! @return true, if active. Otherwise false
  //! @note unit-tested
  inline bool checkIfActive(int pos) const { return m_DofActive.test(pos); }

  //! tell a degree of freedom its position within the system matrices/vectors
  //! @note unit-tested
  inline void setGlobalRow(int dof, PetscInt pos) { m_GlobalRows[dof] = pos; }

  //! read the global position of a degree of freedom
  //! @note unit-tested
  inline PetscInt getGlobalRow(int dof) const { return m_GlobalRows[dof]; }

  //! the whole array - used for copying
  //! @note unit-tested
  inline std::bitset<cstNumberOfKnownDofs> getAllDofFlags(void) const {
    return m_DofActive;
  }

  //! the whole array - used for copying
  //! @note unit-tested
  inline std::vector<PetscInt> getAllGlobalRows(void) const {
    return m_GlobalRows;
  }

  //! reads an object from a stream
  //! @param is inputstream
  //! @return modified inputstream
  std::istream &read(std::istream &is);

  //! writes this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cDegreesOfFreedom &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cDegreesOfFreedom &other) {
  return other.write(os);
}

#endif
