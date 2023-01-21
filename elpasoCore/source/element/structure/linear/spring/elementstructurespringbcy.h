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

#ifndef INFAM_ELEMENTSTRUCTURE_SPRINGBCY_H
#define INFAM_ELEMENTSTRUCTURE_SPRINGBCY_H

#include "../elementstructurelinear.h"

class cElementStructureSpringBCy
    : public cElementStructureLinear,
      private tCounter<cElementStructureSpringBCy> {
 protected:
  //! return values for shape functions evaluated at Gauss points.
  //! Dummy function here.
  cArray3d getShapeFunctions(void) const {
    cArray3d dummy;
    return dummy;
  }

  //! dummy function here
  std::vector<short> getIndicesOfFaceNodes(int) const {
    std::vector<short> dummy(1, -1);
    return dummy;
  }

  //! dummy function here
  short getNumberOfFaces() const { return -1; }

 public:
  cElementStructureSpringBCy();
  cElementStructureSpringBCy(const cElementStructureSpringBCy &other);
  ~cElementStructureSpringBCy();

  //! return the shape of this element
  eElementShape getElementShape(void) const { return Point; }

  //! return the number of faces of this elements (dummy function here)
  short getNumberOfNodesPerFace() const { return -1; }

  //! return the type of element (sometimes better to use this
  //! instead of typeid or dynamic_cast<>
  eElementType getElementType(void) const { return SpringBC; }

  //! return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementStructureSpringBCy>::howMany();
  }

  //! return vector that contains elements degrees of freedom
  std::vector<eKnownDofs> getDofs(void) const;

  //! assembling element's stiffnessmatrix
  //! @param Solution is the calculated global Solution of the System
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector (full vector needed)
  //! @param dx change of x (full vector needed)
  void assembleStiffnessMatrix(cElementMatrix &KM, Vec *x, Vec *dx);

  //! assembling element's massmatrix
  //! @param MM the assembled massmatrix
  void assembleMassMatrix(cElementMatrix &MM);

  //! assembling element's loadvector
  //! @param Solution is the calculated global solution of the system
  //! @param LV the assembled loadvector
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector
  //! @param dx change of x
  void assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x,
                          Vec *dx);

  //! compute stresses element wise
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesCells(Vec &fullSolution, Vec &stresses,
                            Vec &stressesSec2) {
    trace(
        "ERROR: cElementStructureSpringBCy::computeStressesCells() not "
        "implemented, yet.");
    ExitApp();
  }
  //! get stresses element wise without adding to global stresses (needed for
  //! shell element)
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stress2d a vector that holds the local stresses
  //! @param stress2dSec2 a vector that holds the local stresses for section 2
  void getStressesCells(Vec &fullSolution, cElementVector &stresses,
                        cElementVector &stressesSec2) {
    trace(
        "ERROR: cElementStructureSpringBC::computeStressesCells() not "
        "implemented, yet.");
    ExitApp();
  }

  //! compute stresses at the nodes
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses) {
    trace(
        "ERROR: cElementStructureSpringBCy::computeStressesNodes() not "
        "implemented, yet.");
    ExitApp();
  }

  //! returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
    trace(
        "ERROR: "
        "cElementStructureSpringBCy::getLocalCoordinatesOfElementsNodes() not "
        "implemented, yet.");
    ExitApp();
  }

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructureSpringBCy &other) {
  return other.write(os);
}

#endif
