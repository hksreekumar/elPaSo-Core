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

#ifndef INFAM_SHAPE_QUAD4_H
#define INFAM_SHAPE_QUAD4_H

#include "shape.h"

//! @brief class for quadrilaterals with 4 nodes
//! @author Dirk Clasen
//! @date 29.04.2007
class cShapeQuad4 : public cShape, private tCounter<cShapeQuad4> {
 private:
  //! returns the number of class' instances
  static size_t howMany(void) { return tCounter<cShapeQuad4>::howMany(); }

 protected:
  static cArray3d N;  ///< shapefunctions at Gauss points

  static cArray3d N_map;  ///< mapped shapefunctions at Gauss points (only for
                          ///< PoroPlateKienzler)
                          // poro plate KR

  //! evaluated shape functions for this element type at Gauss points
  void initializeShapeFunctions(void);

  void initializeShapeFunctions_N_map(void);  // (only for PoroPlateKienzler) KR

  //! return values of evaluated shapefunctions
  // cArray3d getShapeFunctions(void) const { return N; }

  //! return the number of nodes that describe one face of the element
  short getNumberOfNodesPerFace(void) const { return 2; }

  //! Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..3]
  // std::vector<short> getIndicesOfFaceNodes(int Face) const;

 public:
  cShapeQuad4();
  cShapeQuad4(const cShapeQuad4 &other);
  virtual ~cShapeQuad4();

  //! return values of evaluated shapefunctions
  cArray3d getShapeFunctions(void) const { return N; }

  //! denotes the shape of this element
  eElementShape getElementShape(void) const { return Quadrilateral; }

  //! return the number of faces of this elements
  short getNumberOfFaces(void) const { return 4; }

  //! Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..3]
  std::vector<short> getIndicesOfFaceNodes(int Face) const;
};

#endif
