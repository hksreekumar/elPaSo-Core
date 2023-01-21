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

#ifndef INFAM_MATERIAL_FLUID_IDEAL_H
#define INFAM_MATERIAL_FLUID_IDEAL_H

#include "../materialfluidlin.h"

/**
 * @brief standard linear ideal fluid material
 * @author Dirk Clasen
 */
class cMaterialFluidIdeal :
  public cMaterialFluidLin,
  private tCounter<cMaterialFluidIdeal>
{
public:
  cMaterialFluidIdeal();
  cMaterialFluidIdeal(const cMaterialFluidIdeal &other);
  ~cMaterialFluidIdeal();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cMaterialFluidIdeal>::howMany(); }

  //! compute the frequency dependent wave number.
  //! Its equal to cf for loss free fluids.
  //! @note unit - tested
  PetscScalar getCfOmega(void) const { return getCf(); }

  //! return frequency dependent density
  //! @note unit - tested
  PetscScalar getRhoOmega(void) const { return getRho(); }

  //! update material parameters if they depend on frequency
  void updateMaterial(void) { /* nothing to do */ }

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream& read(std::istream &is);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream& write(std::ostream &os) const;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream& writeXml(std::ostream &os) const;
};


//! overloaded inputoperator
inline std::istream& operator>>(std::istream &is, cMaterialFluidIdeal &other)
{
  return other.read(is);
}


//! overloaded outputoperator
inline std::ostream& operator<<(std::ostream &os, const cMaterialFluidIdeal &other)
{
  return other.write(os);
}

#endif
