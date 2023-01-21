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

#ifndef INFAM_MATERIAL_STRUCTURE_ISO_LIN_H
#define INFAM_MATERIAL_STRUCTURE_ISO_LIN_H

#include "../materialstructureiso.h"

//! @brief linear isotropic material of structural elements
//! @author Marco Schauer
//! @date 04.09.2008
class cMaterialStructureIsoLin : public cMaterialStructureIso,
                                 private tCounter<cMaterialStructureIsoLin> {
 private:
 public:
  cMaterialStructureIsoLin(void);
  cMaterialStructureIsoLin(const cMaterialStructureIsoLin &other);
  ~cMaterialStructureIsoLin();

  //! return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cMaterialStructureIsoLin>::howMany();
  }
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is,
                                cMaterialStructureIsoLin &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cMaterialStructureIsoLin &other) {
  return other.write(os);
}

#endif
