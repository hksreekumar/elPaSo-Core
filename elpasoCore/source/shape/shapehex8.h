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

#ifndef INFAM_SHAPE_HEX8_H
#define INFAM_SHAPE_HEX8_H

#include "shape.h"

//! @brief class for hexahedrons with 8 nodes
//! @author Dirk Clasen
//! @date 29.04.2007
class cShapeHex8 : public cShape, private tCounter<cShapeHex8> {
 private:
  //! returns the number of class' instances
  static size_t howMany(void) { return tCounter<cShapeHex8>::howMany(); }

  static cArray3d Nface;  ///< shapefunctions at Gauss points of element's face

 protected:
  static cArray3d N;  ///< shapefunctions at Gauss points

  //! evaluated shape functions for this element type at Gauss points
  void initializeShapeFunctions(void);

  //! return values of evaluated shapefunctions
  // cArray3d getShapeFunctions(void) const { return N; }

  //! evaluated shapefunctions for element's faces
  cArray3d getShapeFunctionsFace(void) const { return Nface; }

  //! return the number of nodes that describe one face of the element
  short getNumberOfNodesPerFace(void) const { return 4; }

  //! Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..6]
  // std::vector<short> getIndicesOfFaceNodes(int Face) const;

 public:
  cShapeHex8();
  cShapeHex8(const cShapeHex8 &other);
  virtual ~cShapeHex8();

  //! return values of evaluated shapefunctions
  cArray3d getShapeFunctions(void) const { return N; }

  //! denotes the shape of this element
  eElementShape getElementShape(void) const { return Hexahedron; }

  //! return the number of faces of this elements
  short getNumberOfFaces(void) const { return 6; }

  //! Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..6]
  std::vector<short> getIndicesOfFaceNodes(int Face) const;
};

#endif
