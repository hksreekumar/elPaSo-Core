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

#ifndef INFAM_ELEMENT_STRUCTURE_KIRCHHOFF4_H
#define INFAM_ELEMENT_STRUCTURE_KIRCHHOFF4_H

#include "../../../../shape/shapequad4.h"
#include "elementstructurekirchhoff.h"

//! @brief Kirchhoff plate elements with 4 nodes
//! @author Dirk Clasen
//! @date 24.08.2005
//! This element is only suited for rectangular element shapes
class cElementStructureKirchhoff4
    : public cElementStructureKirchhoff,
      public virtual cShapeQuad4,
      private tCounter<cElementStructureKirchhoff4> {
 private:
  //! @brief Numbering of one element's face
  //! @brief The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..3]
  std::vector<short> getIndicesOfFaceNodes(int Face) const {
    return cShapeQuad4::getIndicesOfFaceNodes(Face);
  }

  static cArray3d N;  ///< at Gauss points evaluated shape functions

  //! @brief evaluate shape functions at Gauss points
  void initializeShapeFunctions(void);

  //! @brief return evaluated shape functions
  cArray3d getShapeFunctions(void) const { return N; }

  //! @brief Maybe the element's incidence has to be modified in order to
  //! compute m_Lx as well as m_Ly correctly.
  void rotateElement();

 protected:
  //! @brief return the number of nodes that describe one face of the element
  short getNumberOfNodesPerFace(void) const {
    return cShapeQuad4::getNumberOfNodesPerFace();
  }

 public:
  //! @brief Constructor
  cElementStructureKirchhoff4();
  //! @brief Copy constructor
  cElementStructureKirchhoff4(const cElementStructureKirchhoff4 &other);
  //! @brief Destructor
  ~cElementStructureKirchhoff4();

  //! @brief return the shape of this element
  eElementShape getElementShape(void) const {
    return cShapeQuad4::getElementShape();
  }

  //! @brief return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementStructureKirchhoff4>::howMany();
  }

  //! @brief return the number of faces of this elements
  short getNumberOfFaces(void) const { return cShapeQuad4::getNumberOfFaces(); }

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

  //! @brief returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
    trace("ERROR: getLocalCoordinatesOfElementsNodes() not implemented, yet.");
    ExitApp();
  }

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructureKirchhoff4 &other) {
  return other.write(os);
}

#endif
