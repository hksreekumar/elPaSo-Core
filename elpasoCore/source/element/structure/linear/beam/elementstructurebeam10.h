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

#ifndef INFAM_ELEMENT_STRUCTURE_BEAM_10_H
#define INFAM_ELEMENT_STRUCTURE_BEAM_10_H

#include "elementstructurebeam3d.h"

//! @brief 3D Bernoulli beam element with 2 nodes (no tosion)
//! @author Marco Schauer
//! @date 13.08.2009
//! This element type is based upon cubic test functions - the inner
//! degrees of freedom are eliminated on element level.
//! Theory: KNOTHE, K.; WESSELS, H.: "Finite Elemente", Springer Verlag 1992
//! In constrast to the chapter of the book the rotation of the beam is
//! oriented into the other direction. This notation is used according to
//! Prof. Sandbergs' FE-package CALFEM.
//!
//! The assembly procedure of element's stiffnessmatrix uses a value psi. If
//! the value of psi is equal to 1 the Timoshenko beam becomes a Bernoulli
//! beam.
class cElementStructureBeam10 : public cElementStructureBeam3D,
                                private tCounter<cElementStructureBeam10> {
 private:
  //! @brief computes the transformationmatrix from local to global
  //! coordinate system
  void computeTransformationMatrix(cMatrix &T,
                                   const bool &transpose = false) const;

 protected:
 public:
  //! @brief Constructor
  cElementStructureBeam10(const eUseBeamTheory &Theory = Bernoulli);
  //! @brief Copy constructor
  cElementStructureBeam10(const cElementStructureBeam10 &other);
  //! @brief Destructor
  ~cElementStructureBeam10();

  //! @brief return the number of instances of this object type
  static size_t howMany(void) {
    return tCounter<cElementStructureBeam10>::howMany();
  }

  //! @brief return vector that contains elements degrees of freedom
  std::vector<eKnownDofs> getDofs(void) const;

  //! @brief this function transforms the local normal vector of the beam
  //! into global coordinates. These function is used to define the
  //! coupling direction to acoustic elements
  cVector getGlobalNormalVector(void) const;

  //! @brief assembling element's stiffnessmatrix
  //! @param Solution is the calculated global Solution of the System
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector (full vector needed)
  //! @param dx change of x (full vector needed)
  void assembleStiffnessMatrix(cElementMatrix &KM, Vec *x, Vec *dx);

  //! @brief assembling element's massmatrix
  //! @param MM the assembled massmatrix
  void assembleMassMatrix(cElementMatrix &MM);

  //! @brief assembling element's loadvector
  //! @param Solution is the calculated global solution of the system
  //! @param LV the assembled loadvector
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector
  //! @param dx change of x
  void assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x,
                          Vec *dx);

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief write this object to XML
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructureBeam10 &other) {
  return other.write(os);
}

#endif
