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

#ifndef INFAM_ELEMENT_STRUCTURE_BRICK8_H
#define INFAM_ELEMENT_STRUCTURE_BRICK8_H

#include "../../../../shape/shapehex8.h"
#include "elementstructurebrick.h"

//! @brief solid brick element (8 nodes)
//! @author Dirk Clasen
//! @date 28.08.2005
class cElementStructureBrick8 : public cElementStructureBrick,
                                public virtual cShapeHex8,
                                private tCounter<cElementStructureBrick8> {
 private:
  //! @brief Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..5]
  std::vector<short> getIndicesOfFaceNodes(int Face) const {
    return cShapeHex8::getIndicesOfFaceNodes(Face);
  }

  //! @brief return the array of evaluated shape functions
  cArray3d getShapeFunctions(void) const {
    return cShapeHex8::getShapeFunctions();
  }

  //! @brief return values of evaluated shapefunctions of element's face
  cArray3d getShapeFunctionsFace(void) const {
    return cShapeHex8::getShapeFunctionsFace();
  }

  //! @brief returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz);

 protected:
  //! @brief return the number of nodes that describe one face of the element
  short getNumberOfNodesPerFace(void) const {
    return cShapeHex8::getNumberOfNodesPerFace();
  }

 public:
  //! @brief Constructor
  cElementStructureBrick8();
  //! @brief Copy constructor
  cElementStructureBrick8(const cElementStructureBrick8 &other);
  //! @brief Destructor
  ~cElementStructureBrick8();

  //! @brief return the number of faces of this elements
  short getNumberOfFaces(void) const { return cShapeHex8::getNumberOfFaces(); }

  //! @brief denotes the shape of this element
  eElementShape getElementShape(void) const {
    return cShapeHex8::getElementShape();
  }

  //! @brief return the number of element's instances
  static size_t howMany(void) {
    return tCounter<cElementStructureBrick8>::howMany();
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
                                const cElementStructureBrick8 &other) {
  return other.write(os);
}

#endif
