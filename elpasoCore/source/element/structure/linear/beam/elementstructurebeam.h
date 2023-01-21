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

#ifndef INFAM_ELEMENT_STRUCTURE_BEAM_H
#define INFAM_ELEMENT_STRUCTURE_BEAM_H

#include "../elementstructurelinear.h"

//! @brief Timoshenko/Bernoulli beam element with 2 nodes
//! @author Dirk Clasen
//! @date 28.08.2005
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
class cElementStructureBeam : public cElementStructureLinear,
                              private tCounter<cElementStructureBeam> {
 private:
  eUseBeamTheory m_UseTheory;  ///< specify the beam theory that shall be used

  //! @brief compute the shear parameter psi. If psi==1 we can use this element
  //! formulation to solve Bernoulli beams.
  PetscReal getPsi(void) const;

  //! @brief compute beams length
  inline PetscReal getLength(void) const {
    return m_Nodes[0]->distance(*(m_Nodes[1]));
  }

  //! @brief computes the transformationmatrix from local to global
  //! coordinate system
  void computeTransformationMatrix(cMatrix &T,
                                   const bool &transpose = false) const;

  //! @brief dummy function here
  std::vector<short> getIndicesOfFaceNodes(int) const {
    std::vector<short> dummy(1, -1);
    return dummy;
  }

  //! @brief dummy function here
  short getNumberOfFaces() const { return -1; }

 protected:
  //! @brief return values for shape functions evaluated at Gauss points.
  //! Dummy function here.
  cArray3d getShapeFunctions(void) const {
    cArray3d dummy;
    return dummy;
  }

 public:
  //! @brief Constructor
  cElementStructureBeam(const eUseBeamTheory &Theory = Bernoulli);
  //! @brief Copy constructor
  cElementStructureBeam(const cElementStructureBeam &other);
  //! @brief Destructor
  ~cElementStructureBeam();

  //! @brief return the beam theory used for this element
  eUseBeamTheory getBeamTheory(void) const { return m_UseTheory; }

  //! @brief return the number of faces of this elements (dummy function here)
  short getNumberOfNodesPerFace() const { return -1; }

  //! @brief denotes the shape of this element
  eElementShape getElementShape(void) const { return Beam; }

  //! @brief return the type of element (sometimes better to use this
  //! instead of typeid or dynamic_cast<>
  eElementType getElementType(void) const { return TimoshenkoBeam; }

  //! @brief return the number of instances of this object type
  static size_t howMany(void) {
    return tCounter<cElementStructureBeam>::howMany();
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

  //! @brief compute stresses element wise
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesCells(Vec &fullSolution, Vec &stresses,
                            Vec &stressesSec2) {
    trace(
        "ERROR: cElementStructureBeam::computeStressesCells() not implemented, "
        "yet.");
    ExitApp();
  }
  //! @brief get stresses element wise without adding to global stresses (needed
  //! for shell element)
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stress2d a vector that holds the local stresses
  //! @param stress2dSec2 a vector that holds the local stresses for section 2
  void getStressesCells(Vec &fullSolution, cElementVector &stresses,
                        cElementVector &stressesSec2) {
    trace(
        "ERROR: cElementStructureBeam::getStressesCells() not implemented, "
        "yet.");
    ExitApp();
  }
  //! @brief compute stresses at the nodes
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses) {
    trace(
        "ERROR: cElementStructureBeam::computeStressesNodes() not implemented, "
        "yet.");
    ExitApp();
  }

  //! @brief returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
    trace("ERROR: getLocalCoordinatesOfElementsNodes() not implemented, yet.");
    ExitApp();
  }

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
                                const cElementStructureBeam &other) {
  return other.write(os);
}

#endif
