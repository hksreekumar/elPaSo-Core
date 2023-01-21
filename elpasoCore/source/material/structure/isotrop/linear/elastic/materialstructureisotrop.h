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

#ifndef INFAM_MATERIAL_STRUCTURE_ISOTROP_H
#define INFAM_MATERIAL_STRUCTURE_ISOTROP_H

#include "../materialstructureisolin.h"

//! @brief isotropic material of structural elements
//! @author Dirk Clasen
//! @date 22.06.2005
class cMaterialStructureIsotrop : public cMaterialStructureIsoLin,
                                  private tCounter<cMaterialStructureIsotrop> {
 private:
  PetscScalar m_E;  ///< Young's modulus
  PetscReal m_Nu;   ///< Poisson's ratio

 public:
  cMaterialStructureIsotrop(void);
  cMaterialStructureIsotrop(const cMaterialStructureIsotrop &other);
  ~cMaterialStructureIsotrop();

  //! return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cMaterialStructureIsotrop>::howMany();
  }

  //! return the value of the Young's modulus
  //! @see testmaterialstructureisotrop.h
  PetscScalar getE(void) const { return m_E; }

  //! return the value of the frequency dependend value
  //! of the Young's modulus
  //! @see testmaterialstructureisotrop.h
  PetscScalar getEOmega(void) const { return getE(); }

  //! return Poisson's ratio
  //! @see testmaterialstructureisotrop.h
  PetscReal getNu(void) const { return m_Nu; }

  //! spring constant for structural interface elements
  //! @see testmaterialstructureisotrop.h
  PetscScalar getSpringConst(void) const { return m_E; }

  //! update the material parameters that are frequency dependent
  void updateMaterial(void) {}

  //! return the thickness; new: cElement *elemPtr = nullptr
  //! @see testmaterialstructureisotrop.h
  PetscReal getT(cElement *elemPtr = NULL) const { return m_T; }

  //! shear part of the elasticity matrix of Mindlin plate
  //! @see testmaterialstructureisotrop.h
  void setupCs(cElementMatrix &Cs, cElement *elemPtr = NULL) const;

  //! bending part of the elasticity matrix of Mindlin plate
  //! @see testmaterialstructureisotrop.h
  void setupCb(cElementMatrix &Cb, cElement *elemPtr = NULL) const;

  /*BEGIN_NO_COVERAGE*/
  //! prestress part of the elasticity matrix of Mindlin plate
  void setupCpre(cElementMatrix &Cpre, cElement *elemPtr = NULL) const {}
  /*END_NO_COVERAGE*/

  //! elasticity matrix of membrane elements
  //! @see testmaterialstructureisotrop.h
  void setupCm(cElementMatrix &Cm) const;

  /*BEGIN_NO_COVERAGE*/
  //! elasticity matrix of 3d solid elements
  void setupC(cElementMatrix &C) const;

  //! elasticity matrix for plane stress problem
  void setupCps(cElementMatrix &Cps) const;
  /*END_NO_COVERAGE*/

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
inline std::istream &operator>>(std::istream &is,
                                cMaterialStructureIsotrop &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cMaterialStructureIsotrop &other) {
  return other.write(os);
}

#endif
