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

#ifndef INFAM_ELEMENT_FLUID_H
#define INFAM_ELEMENT_FLUID_H

#include "../../bc/boundaryconditionfluid.h"
#include "../../material/fluid/materialfluid.h"
#include "../elementfem.h"

//! @brief base class for all fluidelements
//! @author Dirk Clasen
//! @date 12.08.2005
class cElementFluid : public cElementFEM, private tCounter<cElementFluid> {
 private:
  //! @brief evaluate impedance boundary condition
  //! This function will simply call subEvaluateImpedance
  //! @param omega angular frequency in s^{-1}
  //! @param KM    dynamic stiffnessmatrix
  virtual void evaluateImpedance(const PetscReal &omega,
                                 cElementMatrix &KM) = 0;

 protected:
  cMaterialFluid *m_Material;  ///< set of materialparameters

  //! @brief return the number of nodes that describe one face of the element
  virtual short getNumberOfNodesPerFace(void) const = 0;

 public:
  //! @brief Constructor
  cElementFluid(short NumberOfNodes, short NumberOfGaussPoints);
  //! @brief Copy constructor
  cElementFluid(const cElementFluid &other);
  //! @brief Destructor
  virtual ~cElementFluid();

  //! @brief returns the physical meaning of this element type
  ePhysicsType getPhysicsType(void) const { return Acoustics; }

  //! @brief return the type of element (sometimes better to use this
  //! instead of typeid or dynamic_cast<>
  eElementType getElementType(void) const { return Helmholtz; }

  //! @brief return the number of faces of this elements
  virtual short getNumberOfFaces(void) const = 0;

  //! @brief return the normal vector on one face of the 2d fluid element.
  //! The normal vector points outward of the element.
  virtual cVector getNormalOnFace(const short &face) const = 0;

  //! @brief returns the number of class' instances
  static size_t howMany(void) { return tCounter<cElementFluid>::howMany(); }

  //! @brief return vector that contains elements degrees of freedom
  std::vector<eKnownDofs> getDofs(void) const;

  //! @brief assemble dynamic stiffness matrix (frequency domain)
  //! In general, K-omega^2 M will be computed. For fluidelements e.g. impedance
  //! boundary conditions will be incorporated.
  //! @param omega angular frequency omega = 2 \pi f [s^{-1}]
  //! @param EM dynamic stiffnessmatrix
  void assembleDynamicStiffnessMatrix(const PetscReal &omega,
                                      cElementMatrix &EM);

  //! @brief assemble frequency dependent loadvector
  //! @param omega  angular frequency omega = 2 \pi f [s^{-1}]
  //! @param LV  dynamic elementloadvector
  virtual void assembleDynamicLoadVector(const PetscReal &omega,
                                         cElementVector &LV,
                                         cElementMatrix &EM) = 0;

  //! @brief compute stresses element wise
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesCells(Vec &fullSolution, Vec &stresses,
                            Vec &stressesSec2) {
    trace("ERROR: cElementFluid::computeStressesCells() not implemented, yet.");
    ExitApp();
  }
  //! @brief get stresses element wise without adding to global stresses (needed
  //! for shell element)
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stress2d a vector that holds the local stresses
  //! @param stress2dSec2 a vector that holds the local stresses for section 2
  void getStressesCells(Vec &fullSolution, cElementVector &stresses,
                        cElementVector &stressesSec2) {
    trace("ERROR: cElementFluid::getStressesCells() not implemented, yet.");
  }
  //! @brief compute stresses at the nodes
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses) {
    trace("ERROR: cElementFluid::computeStressesNodes() not implemented, yet.");
    ExitApp();
  }

  //! @brief returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
    trace("ERROR: getLocalCoordinatesOfElementsNodes() not implemented, yet.");
    ExitApp();
  }

  //! @brief return material properties
  cMaterial *getMaterial(void) const;

  //! @brief Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..5]
  virtual std::vector<short> getIndicesOfFaceNodes(int Face) const = 0;

  //! assign a set of materialparameters to the element
  void setMaterial(cMaterial *ptrMaterial);

  //! @brief writes this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief writes this object in XML
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded output operator
inline std::ostream &operator<<(std::ostream &os, const cElementFluid &other) {
  return other.write(os);
}

#endif
