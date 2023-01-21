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

#ifndef INFAM_ELEMENTLOAD_FLUID_H
#define INFAM_ELEMENTLOAD_FLUID_H

#include "elementload.h"

//! @brief base class for element loads of fluid elements
class cElementLoadFluid : public cElementLoad,
                          private tCounter<cElementLoadFluid> {
 private:
  PetscReal m_Source;

 public:
  //! @brief Constructor
  cElementLoadFluid();
  //! @brief Copy constructor
  cElementLoadFluid(const cElementLoadFluid &other);
  //! @brief Destructor
  ~cElementLoadFluid();

  //! @brief number of instances of this object type
  static size_t howMany(void) { return tCounter<cElementLoadFluid>::howMany(); }

  inline PetscReal getSourceValue(void) const { return m_Source; }

  std::ostream &write(std::ostream &os) const;
  std::ostream &writeXml(std::ostream &os) const;
  std::istream &read(std::istream &is);
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementLoadFluid &other) {
  return other.write(os);
}

//! @brief overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cElementLoadFluid &other) {
  return other.read(is);
}

#endif
