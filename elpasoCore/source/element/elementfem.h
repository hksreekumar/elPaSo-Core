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

#ifndef INFAM_ELEMENT_FEM_H
#define INFAM_ELEMENT_FEM_H

#include "../shape/shape.h"
#include "element.h"

typedef std::multimap<short, cElementLoad *> MapLoadsOnElement;
typedef MapLoadsOnElement::const_iterator ItLoadsOnElementMap;

//! @brief virtual base class of all classes of finite elements
//! @author Dirk Clasen
//! @date 15.09.2005
class cElementFEM : public cElement,
                    public virtual cShape,
                    private tCounter<cElementFEM> {
 private:
  const short m_NumberOfDofsPerNode;  ///< number of degrees of freedom per node
  const short m_NumberOfGaussPoints;  ///< number of Gauss points used for
                                      ///< numerical integration

 protected:
  static cMatrix Jac;  ///< Jacobian
  static cMatrix Jac1;
  static PetscReal detJac;  ///< determinant of the Jacobian
  static cMatrix invJac;    ///< inverse of Jacobian

  MapLoadsOnElement m_ElementLoads;  ///< elementloads

  //! @brief compute the Jacobian as well as its determinant (2D case)
  //! @param N  at all Gauss points evaluated test functions
  //! @param gp number of the Gauss point for which Jac has to be computed
  void setupJacobian2D(cArray3d N, int gp);

  //! @brief compute the Jacobian as well as its determinant (2D Plate case)
  //! @param N  at all Gauss points evaluated test functions
  //! @param gp number of the Gauss point for which Jac has to be computed
  //! only to be used in Kienzler poro plate
  void setupJacobian2D_plate(cArray3d N, int gp, cArray3d N_map);  // KR

  //! @brief compute the Jacobian as well as its determinant (3D case)
  //! @param N  at all Gauss points evaluated test functions
  //! @param gp number of the Gauss point for which Jac has to be computed
  void setupJacobian3D(cArray3d N, int gp);
  // void setupJacobian3DL(cArray3d N, int gp);

  //! @brief invert the Jacobian (2D)
  void invertJacobian2D(void);

  //! @brief invert the Jacobian (3D)
  void invertJacobian3D(void);
  // void invertJacobian3DL(void);

  //! @brief returns the number of instances of the current object
  static size_t howMany(void) { return tCounter<cElementFEM>::howMany(); }

 public:
  std::string m_Orientation;  ///< orientation of material in this element
  std::vector<PetscReal>
      m_OrientationVector;  ///< orientation vector of material in this element,
                            ///< if m_Orientation is "user-def"

  cElementFEM(short NumberOfNodes, short NumberOfDofsPerNode,
              short NumberOfGaussPoints);
  cElementFEM(const cElementFEM &other);
  virtual ~cElementFEM();

  //! @brief denotes the shape of this element
  virtual eElementShape getElementShape(void) const = 0;

  //! @brief return the type of element (sometimes better to use this
  //! instead of typeid or dynamic_cast<>
  virtual eElementType getElementType(void) const = 0;

  //! @brief return the number of degrees of freedom per elementnode
  inline short getNumberOfDofsPerNode(void) const {
    return m_NumberOfDofsPerNode;
  }

  //! @brief return the number of Gauss points to use for numerical integration
  inline short getNumberOfGaussPoints(void) const {
    return m_NumberOfGaussPoints;
  }

  //! @brief return vector that contains elements degrees of freedom
  virtual std::vector<eKnownDofs> getDofs(void) const = 0;

  //! @brief Returns global numbering of element's degrees of freedom.
  //! Used for assembling global matrices and vectors
  std::vector<PetscInt> getGlobalPositions(void) const;

  //! apply nodal boundary conditions to element matrix
  // void applyNodalBoundaryConditions(cElementMatrix &KM, cElementVector &LV);

  //! @brief assign an elementload to this element
  //! @param face surface to which the load has to be applied (only important
  //! for volume elements)
  //! @param ptrElementLoad points to an elementload
  void insertElementLoad(const int &face, cElementLoad *ptrElementLoad);

  //! @brief assembling element's stiffnessmatrix
  //! @param Solution is the calculated global Solution of the System
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector (full vector needed)
  //! @param dx change of x (full vector needed)
  virtual void assembleStiffnessMatrix(cElementMatrix &KM, Vec *x, Vec *dx) = 0;

  //! @brief assembling element's massmatrix
  //! @param MM the assembled massmatrix
  virtual void assembleMassMatrix(cElementMatrix &MM) = 0;

  //! @brief assembling element's loadvector
  //! @param Solution is the calculated global solution of the system
  //! @param LV the assembled loadvector
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector
  //! @param dx change of x
  virtual void assembleLoadVector(cElementVector &LV, cElementMatrix &KM,
                                  Vec *x, Vec *dx) = 0;

  //! @brief extract local solution from full solution
  //! @param fullSolution is the calculated global solution of the system
  //! @param localSolution the assembled local solution
  void getLocalSolutionFromFullSolution(Vec *fullSolution,
                                        cElementVector &localSolution);

  //! @brief computes the dynamic stiffness matrix for calculations in
  //! frequency domain. In general, only K - omega^2 M will be computed.
  //! For fluid elements, some boundary conditions (impedance) will be applied.
  //! @param omega angular frequency omega = 2 \pi f [s^{-1}]
  //! @param EM dynamic stiffnessmatrix
  virtual void assembleDynamicStiffnessMatrix(const PetscReal &omega,
                                              cElementMatrix &EM) = 0;

  //! @brief computes dynamic load vector
  //! @param omega  angular frequency omega = 2 \pi f [s^{-1}]
  //! @param LV  dynamic loadvector of this element
  //! @param EM dynamic stiffnessmatrix
  virtual void assembleDynamicLoadVector(const PetscReal &omega,
                                         cElementVector &LV,
                                         cElementMatrix &EM) = 0;

  //! @brief compute stresses element wise
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stresses a vector that holds the global stresses
  virtual void computeStressesCells(Vec &fullSolution, Vec &stresses,
                                    Vec &stressesSec2) = 0;

  //! @brief get stresses element wise without adding to global stresses (needed
  //! for shell element)
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stress2d a vector that holds the local stresses
  //! @param stress2dSec2 a vector that holds the local stresses for section 2
  // virtual void getStressesCells(Vec &fullSolution, cElementVector &stress2d,
  // cElementVector &stress2dSec2) = 0;

  //! @brief compute stresses at the nodes
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stresses a vector that holds the global stresses
  virtual void computeStressesNodes(Vec &fullSolution, Vec &stresses) = 0;

  //! @brief returns local coordinates of element's nodes
  virtual void getLocalCoordinatesOfElementsNodes(cMatrix &xyz) = 0;

  //! @brief returns material properties of this element
  virtual cMaterial *getMaterial(void) const = 0;

  //! @brief assings material properties to this element.
  //! It is also checked if the set of material parameters is
  //! appropriate for the current element.
  //! @param ptrMaterial points to a set of material parameters
  virtual void setMaterial(cMaterial *ptrMaterial) = 0;

  //! @brief sets orientation type for this element (string)
  void setOrientation(std::string orientation) { m_Orientation = orientation; }

  //! @brief sets orientation vector for this element (vector, only available if
  //! orientation type is "user-def")
  void setOrientationVector(std::vector<PetscReal> orientation) {
    m_OrientationVector = orientation;
  }

  //! @brief returns the number of element loads assigned
  short getNumberOfLoads(void) const { return (short)m_ElementLoads.size(); }

  //! @brief returns an iterator to the first elementload of this element
  ItLoadsOnElementMap getFirstElementLoad(void) const {
    return m_ElementLoads.begin();
  }

  //! @brief returns an iterator to the last elementload of this element
  ItLoadsOnElementMap getLastElementLoad(void) const {
    return m_ElementLoads.end();
  }

  //! @brief return the ids of the nodes that describe a face of the element
  std::vector<PetscInt> getIdsOfFaceNodes(int Face) const;

  //! @brief calculates internal forces and deducts them from loadvector
  void assembleInternalForces(Vec *fullSolution, cElementMatrix &KM,
                              cElementVector &LV);

  //! @brief Writes this object to a stream.
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief Writes this object in XML format.
  //! This function will be called after looking for interface elements.
  //! @param os outputstream
  //! @return der modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cElementFEM &other) {
  return other.write(os);
}

#endif
