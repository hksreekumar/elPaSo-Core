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

#ifndef INFAM_ELEMENT_STRUCTURE_MINDLIN_DSG4_H
#define INFAM_ELEMENT_STRUCTURE_MINDLIN_DSG4_H

#include "../../../../shape/shapequad4.h"
#include "elementstructuremindlin.h"

//! @brief Mindlin plate element with 4 nodes
//! @author Dirk Clasen
//! @date 09.08.2005
//! quadrilateral 4 node plate element with bilinear test functions.
//! In order to eliminate the so-called shear locking the DSG
//! approach is eliminated.
class cElementStructureMindlinDSG4
    : public cElementStructureMindlin,
      public virtual cShapeQuad4,
      private tCounter<cElementStructureMindlinDSG4> {
 private:
  //! @brief Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..3]
  std::vector<short> getIndicesOfFaceNodes(int Face) const {
    return cShapeQuad4::getIndicesOfFaceNodes(Face);
  }

  //! @brief return array with the values for the shape functions
  //! computed at the Gauss points
  cArray3d getShapeFunctions(void) const {
    return cShapeQuad4::getShapeFunctions();
  }

 protected:
  //! @brief return the number of nodes that describe one face of the element
  short getNumberOfNodesPerFace(void) const {
    return cShapeQuad4::getNumberOfNodesPerFace();
  }

 public:
  //! @brief Constructor
  cElementStructureMindlinDSG4(void);
  //! @brief Copy constructor
  cElementStructureMindlinDSG4(const cElementStructureMindlinDSG4 &other);
  //! @brief Destructor
  ~cElementStructureMindlinDSG4();

  //! @brief return the shape of this element
  eElementShape getElementShape(void) const {
    return cShapeQuad4::getElementShape();
  }

  //! @brief return the number of faces of this elements
  short getNumberOfFaces(void) const { return cShapeQuad4::getNumberOfFaces(); }

  //! @brief return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementStructureMindlinDSG4>::howMany();
  }

  //! @brief assembling element's stiffnessmatrix
  //! @param Solution is the calculated global Solution of the System
  //! @param KM the assembled stiffnessmatrix
  //! @param x solution vector (full vector needed)
  //! @param dx change of x (full vector needed)
  void assembleStiffnessMatrix(cElementMatrix &KM, Vec *x, Vec *dx);

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
                                const cElementStructureMindlinDSG4 &other) {
  return other.write(os);
}

#endif
