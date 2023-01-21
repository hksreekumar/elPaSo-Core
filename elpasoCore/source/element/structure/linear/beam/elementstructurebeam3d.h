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

#ifndef INFAM_ELEMENT_STRUCTURE_BEAM3D_H
#define INFAM_ELEMENT_STRUCTURE_BEAM3D_H

#include "../elementstructurelinear.h"

//! @brief super class 3D beam element with 2 nodes
//! @author Marco Schauer
//! @date 17.08.2009
class cElementStructureBeam3D : public cElementStructureLinear,
                                private tCounter<cElementStructureBeam3D> {
 private:
  //! @brief computes the transformationmatrix from local to global
  //! coordinate system
  virtual void computeTransformationMatrix(
      cMatrix &T, const bool &transpose = false) const = 0;

 protected:
  //! @brief dummy function here
  std::vector<short> getIndicesOfFaceNodes(int) const {
    std::vector<short> dummy(1, -1);
    return dummy;
  }

  //! @brief dummy function here
  short getNumberOfFaces() const { return -1; }

  eUseBeamTheory m_UseTheory;  ///< specify the beam theory that shall be used

  //! @brief compute beams length
  inline PetscReal getLength(void) const {
    return m_Nodes[0]->distance(*(m_Nodes[1]));
  }

  //! @brief compute the shear parameter psi. If psi==1 we can use this element
  //! formulation to solve Bernoulli beams.
  PetscReal getPsi(void) const;

  //! @brief return values for shape functions evaluated at Gauss points.
  //! Dummy function here.
  cArray3d getShapeFunctions(void) const {
    cArray3d dummy;
    return dummy;
  }

 public:
  //! @brief Constructor
  cElementStructureBeam3D(short NumberOfNodes, short NumberOfDofsPerNode,
                          short NumberOfGaussPoints,
                          const eUseBeamTheory &Theory = Bernoulli);
  //! @brief Copy constructor
  cElementStructureBeam3D(const cElementStructureBeam3D &other);
  //! @brief Destructor
  ~cElementStructureBeam3D();

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
    return tCounter<cElementStructureBeam3D>::howMany();
  }

  //! @brief this function transforms the local normal vector of the beam
  //! into global coordinates. These function is used to define the
  //! coupling direction to acoustic elements
  cVector getGlobalNormalVector(void) const;

  //! @brief compute stresses element wise
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesCells(Vec &fullSolution, Vec &stresses,
                            Vec &stressesSec2) {
    trace(
        "ERROR: cElementStructureBeam3D::computeStressesCells() not "
        "implemented, yet.");
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
        "ERROR: cElementStructureBeam3D::getStressesCells() not implemented, "
        "yet.");
  }
  //! @brief compute stresses at the nodes
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses) {
    trace(
        "ERROR: cElementStructureBeam3D::computeStressesNodes() not "
        "implemented, yet.");
    ExitApp();
  }

  //! @brief returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
    trace("ERROR: getLocalCoordinatesOfElementsNodes() not implemented, yet.");
    ExitApp();
  }
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructureBeam3D &other) {
  return other.write(os);
}

#endif
