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

#ifndef INFAM_ELEMENT_STRUCTURE_LINEAR_H
#define INFAM_ELEMENT_STRUCTURE_LINEAR_H

#include "../../../bc/boundaryconditionstructure.h"
#include "../../../material/structure/materialstructure.h"
#include "../../load/elementloadstructure.h"
#include "../elementstructure.h"

//! @brief virtual base class for linear structural elements (non-poro)
//! @author Dirk Clasen
//! @date 22.05.2006
class cElementStructureLinear : public cElementStructure,
                                private tCounter<cElementStructureLinear> {
 protected:
  cMaterialStructure *m_Material;  ///< material parameters

 public:
  cElementStructureLinear(short NumberOfNodes, short NumberOfDofsPerNode,
                          short NumberOfGaussPoints);
  cElementStructureLinear(const cElementStructureLinear &other);
  virtual ~cElementStructureLinear();

  //! return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cElementStructureLinear>::howMany();
  }

  //! return the set of material parameters assigned to this element
  cMaterial *getMaterial(void) const { return m_Material; }

  //! assign a set of materialparameters appropriate to this element type
  void setMaterial(cMaterial *ptrMaterial);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;
};

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementStructureLinear &other) {
  return other.write(os);
}

#endif
