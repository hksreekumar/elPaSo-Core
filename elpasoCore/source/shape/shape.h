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

#ifndef INFAM_SHAPE_H
#define INFAM_SHAPE_H

#include "../element/element.h"
#include "../misc/array3d.h"
#include "../misc/counter.h"
#include "../misc/mytypes.h"

//! @brief virtual base class for all element shapes
//! @author Dirk Clasen
//! @date 29.04.2007
//! the aim of this approach is to reduce the occurrence of
//! duplicate code used to compute shape functions
class cShape : public virtual cLogging {
 protected:
  //! return values of evaluated shapefunctions
  // virtual cArray3d getShapeFunctions(void) const = 0;

  //! return the number of nodes that describe one face of the element
  virtual short getNumberOfNodesPerFace(void) const = 0;

  //! Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..3]
  // virtual std::vector<short> getIndicesOfFaceNodes(int Face) const = 0;

 public:
  cShape();
  cShape(const cShape &other);
  virtual ~cShape();

  //! return values of evaluated shapefunctions
  virtual cArray3d getShapeFunctions(void) const = 0;

  //! denotes the shape of this element
  virtual eElementShape getElementShape(void) const = 0;

  //! return the number of faces of this elements
  virtual short getNumberOfFaces(void) const = 0;

  //! Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..3]
  virtual std::vector<short> getIndicesOfFaceNodes(int Face) const = 0;
};

#endif
