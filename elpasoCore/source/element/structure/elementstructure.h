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

#ifndef INFAM_ELEMENT_STRUCTURE_H
#define INFAM_ELEMENT_STRUCTURE_H

#include "../elementfem.h"

//! @brief enum for representing the used beam theory
enum eUseBeamTheory { Bernoulli = 1, Timoshenko = 2 };

//! @brief virtual base class for structural elements
//! @author Dirk Clasen
//! @date 22.06.2005
class cElementStructure : public cElementFEM,
                          private tCounter<cElementStructure> {
 public:
  cElementStructure(short NumberOfNodes, short NumberOfDofsPerNode,
                    short NumberOfGaussPoints);
  cElementStructure(const cElementStructure &other);
  virtual ~cElementStructure();

  //! returns the physical meaning of this element type
  ePhysicsType getPhysicsType(void) const { return Elastodynamics; }

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cElementStructure>::howMany(); }

  //! assembling element's dynamic stiffnessmatrix (frequency domain)
  //! In general, EM = (K-omega^2*M) is computed
  //! @param omega angular frequency omega = 2 \pi f [s^{-1}]
  //! @param EM element's dynamic stiffnessmatrix
  void assembleDynamicStiffnessMatrix(const PetscReal &omega,
                                      cElementMatrix &EM);

  //! assemble frequency dependent loadvector
  //! @param omega  angular frequency omega = 2 \pi f [s^{-1}]
  //! @param LV dynamic load vector
  //! @param EM element's dynamic stiffnessmatrix
  void assembleDynamicLoadVector(const PetscReal &omega, cElementVector &LV,
                                 cElementMatrix &EM);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;
};

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructure &other) {
  return other.write(os);
}

#endif
