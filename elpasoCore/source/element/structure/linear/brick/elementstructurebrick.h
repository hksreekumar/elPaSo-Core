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

#ifndef INFAM_ELEMENT_STRUCTURE_BRICK
#define INFAM_ELEMENT_STRUCTURE_BRICK

#include "../elementstructurelinear.h"

//! @brief virtual base for solid elements
//! @author Dirk Clasen
//! @date 28.08.2005
class cElementStructureBrick : public cElementStructureLinear,
                               private tCounter<cElementStructureBrick> {
 private:
  //! @brief computes the surface integral on a face of a hexaedron. The result
  //! is obtained in global coordinates and has to be transformed to local
  //! coordinate system
  //! @param Face surface on which the integral will be computed. valid range:
  //! [0..5]
  //! @param Nface shape functions of element's face
  //! @param C results of the integration
  void evaluateSurfaceIntegral(int Face, cArray3d &Nface, cMatrix &C);

 protected:
  //! @brief evaluated shapefunctions for element's faces
  virtual cArray3d getShapeFunctionsFace(void) const = 0;

  //! @brief return the number of nodes that describe one face of the element
  virtual short getNumberOfNodesPerFace(void) const = 0;

 public:
  //! @brief Constructor
  cElementStructureBrick(short NumberOfNodes, short NumberOfGaussPoints);
  //! @brief Copy constructor
  cElementStructureBrick(const cElementStructureBrick &other);
  //! @brief Destructor
  ~cElementStructureBrick();

  //! @brief return the type of element (sometimes better to use this
  //! instead of typeid or dynamic_cast<>
  eElementType getElementType(void) const { return Continuum; }

  //! @brief return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementStructureBrick>::howMany();
  }

  //! @brief return vector that contains elements degrees of freedom
  std::vector<eKnownDofs> getDofs(void) const;

  //! @brief compute element's stiffnessmatrix
  //! @param KM element's stiffnessmatrix
  //! @param x solution vector (full vector needed)
  //! @param dx change of x (full vector needed)
  void assembleStiffnessMatrix(cElementMatrix &KM, Vec *x, Vec *dx);

  //! @brief compute element's massmatrix
  //! @param KM element's massmatrix
  void assembleMassMatrix(cElementMatrix &MM);

  //! @brief compute element's loadvector
  //! @param LV element's loadvector
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector
  //! @param dx change of x
  void assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x,
                          Vec *dx);

  //! @brief compute stresses element wise
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesCells(Vec &fullSolution, Vec &stresses,
                            Vec &stressesSec2);
  //! @brief get stresses element wise without adding to global stresses (needed
  //! for shell element)
  //! @param fullSolution the full solution vector, each process holds a copy
  //! @param stress2d a vector that holds the local stresses
  //! @param stress2dSec2 a vector that holds the local stresses for section 2
  void getStressesCells(Vec &fullSolution, cElementVector &stresses,
                        cElementVector &stressesSec2) {
    trace(
        "ERROR: cElementStructureBrick::getStressesCells() not implemented, "
        "yet.");
    ExitApp();
  }

  //! @brief compute stresses at the nodes
  //! @param FullSolution the full solution vector, each process holds a copy
  //! @param Stresses a vector that holds the global stresses
  void computeStressesNodes(Vec &fullSolution, Vec &stresses);

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructureBrick &other) {
  return other.write(os);
}

#endif
