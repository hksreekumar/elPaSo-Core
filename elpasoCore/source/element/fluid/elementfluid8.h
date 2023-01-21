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

#ifndef INFAM_ELEMENT_FLUID8_H
#define INFAM_ELEMENT_FLUID8_H

#include "../../shape/shapehex8.h"
#include "elementfluid3d.h"

//! @brief fluidelement with 8 nodes (hexahedron)
//! @author Dirk Clasen
//! @date 12.08.2005
class cElementFluid8 : public cElementFluid3d,
                       public virtual cShapeHex8,
                       private tCounter<cElementFluid8> {
 private:
  //! @brief Numbering of one element's face
  //! The normal vector points outward.
  //! @param Face number of the face. Valid range: [0..5]
  std::vector<short> getIndicesOfFaceNodes(int Face) const {
    return cShapeHex8::getIndicesOfFaceNodes(Face);
  }

  //! @brief returns local coordinates of element's nodes
  void getLocalCoordinatesOfElementsNodes(cMatrix &xyz);

  //! @brief returns the number of class' instances
  static size_t howMany(void) { return tCounter<cElementFluid8>::howMany(); }

  //! @brief return values of evaluated shapefunctions
  cArray3d getShapeFunctions(void) const {
    return cShapeHex8::getShapeFunctions();
  }

  //! @brief return values of evaluated shapefunctions of element's face
  cArray3d getShapeFunctionsFace(void) const {
    return cShapeHex8::getShapeFunctionsFace();
  }

 protected:
  //! @brief return the number of nodes that describe one face of the element
  short getNumberOfNodesPerFace(void) const {
    return cShapeHex8::getNumberOfNodesPerFace();
  }

 public:
  //! @brief Constructor
  cElementFluid8();
  //! @brief Copy constructor
  cElementFluid8(const cElementFluid8 &other);
  //! @brief Destructor
  ~cElementFluid8();

  //! @brief denotes the shape of this element
  eElementShape getElementShape(void) const {
    return cShapeHex8::getElementShape();
  }

  //! @brief return the number of faces of this elements
  short getNumberOfFaces(void) const { return cShapeHex8::getNumberOfFaces(); }

  void computeVi(PetscReal &nimp, PetscScalar &Zn, const Vec &FullSolution,
                 const Vec &Stresses);

  //! @brief writes this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief writes this object in XML
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! @brief overloaded output operator
inline std::ostream &operator<<(std::ostream &os, const cElementFluid8 &other) {
  return other.write(os);
}

#endif
