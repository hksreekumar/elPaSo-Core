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

#ifndef INFAM_MATERIAL_FLUID_H
#define INFAM_MATERIAL_FLUID_H

#include "../material.h"

//! @brief virtual baseclass of all fluid materials
//! @author Dirk Clasen
//! @date 22.06.2005
class cMaterialFluid : public cMaterial, private tCounter<cMaterialFluid> {
 protected:
  PetscReal m_Cf;  ///< speed of sound
  PetscReal m_T;   ///< thickness of the fluid domain (for 2d computations)

 public:
  cMaterialFluid();
  cMaterialFluid(const cMaterialFluid &other);
  ~cMaterialFluid();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cMaterialFluid>::howMany(); }

  //! compute the frequency dependent wave number.
  //! Its equal to getCf() for loss free fluids.
  virtual PetscScalar getCfOmega(void) const = 0;

  //! compute frequency dependent fluid density
  virtual PetscScalar getRhoOmega(void) const = 0;

  //! returns the speed of sound
  inline PetscReal getCf() const { return m_Cf; }

  //! return the thickness of a 2d fluid
  //! @note unit - tested
  inline PetscReal getT(void) const { return m_T; }

  //! assign a new thickness to the fluid material
  //! @note unit - tested
  void setT(const PetscReal &T) { m_T = T; }

  //! assign a new speed of sound to the material
  //! @param newCf
  inline void setCf(const PetscReal &newCf) { m_Cf = newCf; }

  /*BEGIN_NO_COVERAGE*/
  //! return rho_x, rho_y and rho_z values for a given point with
  //! x,y,z-coordinates
  void computeRhoValues(cPoint &point, cMatrix &rho_comp){/* do nothing */};
  /*END_NO_COVERAGE*/

  //! return rho_c for a given point with x,y,z-coordinates
  //! and returns the bulkmodulus
  //! @note unit - tested
  PetscReal getRhoLambda(cPoint &point, PetscReal &rho_c) {
    return 0; /* do nothing */
  };

  /*BEGIN_NO_COVERAGE*/
  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cMaterialFluid &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cMaterialFluid &other) {
  return other.write(os);
}

#endif
/*END_NO_COVERAGE*/