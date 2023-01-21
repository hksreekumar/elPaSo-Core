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

#ifndef INFAM_MATERIAL_STRUCTURE_H
#define INFAM_MATERIAL_STRUCTURE_H

#include "../../element/element.h"
#include "../material.h"

//! @brief virtual base class for the material of structural elements
//! @author Dirk Clasen
//! @date 22.06.2005
class cMaterialStructure : public cMaterial,
                           private tCounter<cMaterialStructure> {
 protected:
  PetscReal m_T;   ///< thickness
  PetscReal m_A;   ///< area of the crosssection
  PetscReal m_Ix;  ///< moment of inertia
  PetscReal m_Iy;
  PetscReal m_Iz;
  PetscReal m_Ks;  ///< shear correction factor
  PetscReal m_Fi;  ///< initial force to prestress element

 public:
  cMaterialStructure();
  cMaterialStructure(const cMaterialStructure &other);
  virtual ~cMaterialStructure();

  //! return value of object counter
  static size_t howMany(void) {
    return tCounter<cMaterialStructure>::howMany();
  }

  //! Young's modulus
  virtual PetscScalar getE(void) const = 0;

  //! return frequency dependent Young's modulus
  virtual PetscScalar getEOmega(void) const = 0;

  //! Poisson's ratio
  virtual PetscReal getNu(void) const = 0;

  //! spring constant for structural interface elements
  virtual PetscScalar getSpringConst(void) const = 0;

  //! spring constant for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscReal getMass(void) const { return -1; }

  //! spring constant for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscScalar getCx(void) const { return -1; }

  //! spring constant for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscScalar getCy(void) const { return -1; }

  //! spring constant for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscScalar getCz(void) const { return -1; }

  //! spring constant for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscScalar getCrx(void) const { return -1; }

  //! spring constant for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscScalar getCry(void) const { return -1; }

  //! spring constant for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscScalar getCrz(void) const { return -1; }

  //! spring lossfactor for structural interface elements
  //! @see testmaterialstructure.h
  virtual PetscReal getEtaSpring(void) const { return -1; }

  //! return crosssection's area
  inline PetscReal getA(void) const { return m_A; }

  //! return the moment of inertia
  //! @see testmaterialstructure.h
  inline PetscReal getI(void) const { return m_Iy; }

  //! return the moment of inertia
  //! @see testmaterialstructure.h
  inline PetscReal getIx(void) const { return m_Ix; }

  //! return the moment of inertia
  //! @see testmaterialstructure.h
  inline PetscReal getIy(void) const { return m_Iy; }

  //! return the moment of inertia
  //! @see testmaterialstructure.h
  inline PetscReal getIz(void) const { return m_Iz; }

  //! return the shear correction factor
  inline PetscReal getKs(void) const { return m_Ks; }

  //! return the thickness; new: cElement *elemPtr = nullptr
  virtual PetscReal getT(cElement *elemPtr = NULL) const = 0;
  // previous implementation: inline PetscReal getT(void) const { return m_T;}

  //! return the thickness
  //! @see testmaterialstructure.h
  inline PetscReal getFi(void) const { return m_Fi; }

  //! assign a new thickness
  void setT(PetscReal T) { m_T = T; }

  //! bending part of the elasticity matrix of Mindlin plate
  virtual void setupCb(cElementMatrix &Cb, cElement *elemPtr = NULL) const = 0;

  //! prestress part of the elasticity matrix of Mindlin plate
  virtual void setupCpre(cElementMatrix &Cpre,
                         cElement *elemPtr = NULL) const = 0;

  //! shear part of the elasticity matrix of Mindlin plate
  virtual void setupCs(cElementMatrix &Cs, cElement *elemPtr = NULL) const = 0;

  //! elasticity matrix of membrane elements
  virtual void setupCm(cElementMatrix &Cm) const = 0;

  //! elasticity matrix of 3d solid elements
  virtual void setupC(cElementMatrix &C) const = 0;

  //! elasticity matrix for plane stress problem
  virtual void setupCps(cElementMatrix &Cps) const = 0;

  /*BEGIN_NO_COVERAGE*/
  //! return rho_x, rho_y and rho_z values for a given point with
  //! x,y,z-coordinates
  void computeRhoValues(cPoint &point, cMatrix &rho_comp){/* do nothing */};
  /*END_NO_COVERAGE*/

  //! return rho_c for a given point with x,y,z-coordinates
  //! and returns the bulkmodulus
  //! @see testmaterialstructure.h
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
inline std::istream &operator>>(std::istream &is, cMaterialStructure &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cMaterialStructure &other) {
  return other.write(os);
}

#endif
/*END_NO_COVERAGE*/