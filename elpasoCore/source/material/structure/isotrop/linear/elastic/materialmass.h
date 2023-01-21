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

#ifndef INFAM_MATERIAL_MASS_H
#define INFAM_MATERIAL_MASS_H

#include "../materialstructureisolin.h"

//! @brief isotropic material of structural elements
//! @author Blech, Rothe
//! @date 24.12.2015
class cMaterialMass : public cMaterialStructureIsoLin,
                      private tCounter<cMaterialMass> {
 private:
  PetscReal m_M;  ///< pointmass

 public:
  cMaterialMass(void);
  cMaterialMass(const cMaterialMass &other);
  ~cMaterialMass();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cMaterialMass>::howMany(); }

  PetscReal getMass(void) const { return m_M; }

  //! return the value of the Young's modulus
  PetscScalar getE(void) const { return 0; }

  //! return the value of the frequency dependend value
  //! of the Young's modulus
  PetscScalar getEOmega(void) const { return getE(); }

  //! return Poisson's ratio
  PetscReal getNu(void) const { return 0; }

  //! spring constant for structural interface elements
  PetscScalar getSpringConst(void) const { return 0; }

  //! update the material paramters that are frequency dependent
  void updateMaterial(void) {}

  //! return the thickness; new: cElement *elemPtr = nullptr
  PetscReal getT(cElement *elemPtr = NULL) const {
    trace(
        "getT() function not defined in isotropic linear elastic mass "
        "material");
    ExitApp();
    return 0.0;
  }

  //! shear part of the elasticity matrix of Mindlin plate
  void setupCs(cElementMatrix &Cs, cElement *elemPtr = NULL) const {}

  //! bending part of the elasticity matrix of Mindlin plate
  void setupCb(cElementMatrix &Cb, cElement *elemPtr = NULL) const {}

  //! prestress part of the elasticity matrix of Mindlin plate
  void setupCpre(cElementMatrix &Cpre, cElement *elemPtr = NULL) const {}

  //! elasticity matrix of membrane elements
  void setupCm(cElementMatrix &Cm) const {}

  //! elasticity matrix of 3d solid elements
  void setupC(cElementMatrix &C) const {}

  //! elasticity matrix for plane stress problem
  void setupCps(cElementMatrix &Cps) const {}

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cMaterialMass &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cMaterialMass &other) {
  return other.write(os);
}

#endif
