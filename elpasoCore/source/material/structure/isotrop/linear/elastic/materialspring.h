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

#ifndef INFAM_MATERIAL_SPRING_H
#define INFAM_MATERIAL_SPRING_H

#include "../materialstructureisolin.h"

//! @brief isotropic material of structural elements
//! @author Blech, Rothe
//! @date 24.12.2015
class cMaterialSpring : public cMaterialStructureIsoLin,
                        private tCounter<cMaterialSpring> {
 private:
  PetscScalar m_Cx;   ///< stiffness trans x
  PetscScalar m_Cy;   ///< stiffness trans y
  PetscScalar m_Cz;   ///< stiffness trans z
  PetscScalar m_Crx;  ///< stiffness rot x
  PetscScalar m_Cry;  ///< stiffness rot y
  PetscScalar m_Crz;  ///< stiffness rot z

  PetscScalar m_Cx_vis;   ///< stiffness trans x
  PetscScalar m_Cy_vis;   ///< stiffness trans y
  PetscScalar m_Cz_vis;   ///< stiffness trans z
  PetscScalar m_Crx_vis;  ///< stiffness rot x
  PetscScalar m_Cry_vis;  ///< stiffness rot y
  PetscScalar m_Crz_vis;  ///< stiffness rot z

  PetscReal m_Eta_Spring;  ///< lossfactor for imaginary part
  short m_ViscoType;       ///< flag that defines the type of damping
#ifdef PETSC_USE_COMPLEX
  std::complex<PetscReal> m_ImagOne;  ///< complex number i=\sqrt{1}
#else

#endif
  std::string m_Filename_Cx;  ///< name of the inputfile for Cx
  PetscInt m_lines_in_Cx_file;
  std::vector<PetscReal> m_Omega_Cx;  ///< frequencies for Cx
  std::vector<PetscReal> m_Cx_var;    ///< Cx, frequency-dependent

  std::string m_Filename_Cy;  ///< name of the inputfile for Cy
  PetscInt m_lines_in_Cy_file;
  std::vector<PetscReal> m_Omega_Cy;  ///< frequencies for Cy
  std::vector<PetscReal> m_Cy_var;    ///< Cy, frequency-dependent

  std::string m_Filename_Cz;  ///< name of the inputfile for Cz
  PetscInt m_lines_in_Cz_file;
  std::vector<PetscReal> m_Omega_Cz;  ///< frequencies for Cz
  std::vector<PetscReal> m_Cz_var;    ///< Cz, frequency-dependent

  std::string m_Filename_Crx;  ///< name of the inputfile for Crx
  PetscInt m_lines_in_Crx_file;
  std::vector<PetscReal> m_Omega_Crx;  ///< frequencies for Crx
  std::vector<PetscReal> m_Crx_var;    ///< Crx, frequency-dependent

  std::string m_Filename_Cry;  ///< name of the inputfile for Cry
  PetscInt m_lines_in_Cry_file;
  std::vector<PetscReal> m_Omega_Cry;  ///< frequencies for Cry
  std::vector<PetscReal> m_Cry_var;    ///< Cry, frequency-dependent

  std::string m_Filename_Crz;  ///< name of the inputfile for Crz
  PetscInt m_lines_in_Crz_file;
  std::vector<PetscReal> m_Omega_Crz;  ///< frequencies for Crz
  std::vector<PetscReal> m_Crz_var;    ///< Crz, frequency-dependent

  std::string m_Filename_Eta_Spring;  ///< name of the inputfile for eta
  PetscInt m_lines_in_Eta_Spring_file;
  std::vector<PetscReal> m_Omega_Eta_Spring;  ///< frequencies for eta
  std::vector<PetscReal> m_Eta_Spring_var;    ///< eta, frequency-dependent

  //! read file storing the values
  void readCxFromFile(void);
  void readCyFromFile(void);
  void readCzFromFile(void);
  void readCrxFromFile(void);
  void readCryFromFile(void);
  void readCrzFromFile(void);
  void readEtaSpringFromFile(void);

 public:
  cMaterialSpring(void);
  cMaterialSpring(const cMaterialSpring &other);
  ~cMaterialSpring();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cMaterialSpring>::howMany(); }

  //! @see testmaterialspring.h
  PetscScalar getCx(void) const { return m_Cx_vis; }

  //! @see testmaterialspring.h
  PetscScalar getCy(void) const { return m_Cy_vis; }

  //! @see testmaterialspring.h
  PetscScalar getCz(void) const { return m_Cz_vis; }

  //! @see testmaterialspring.h
  PetscScalar getCrx(void) const { return m_Crx_vis; }

  //! @see testmaterialspring.h
  PetscScalar getCry(void) const { return m_Cry_vis; }

  //! @see testmaterialspring.h
  PetscScalar getCrz(void) const { return m_Crz_vis; }

  //! @see testmaterialspring.h
  PetscReal getEtaSpring(void) const { return m_Eta_Spring; }

  //! @see testmaterialspring.h
  inline short getType() const { return m_ViscoType; }

  //! @see testmaterialspring.h
  inline void setType(const short &value) { m_ViscoType = value; }

  //! return the value of the Young's modulus
  //! @see testmaterialspring.h
  PetscScalar getE(void) const { return 0; }

  //! return the value of the frequency dependend value
  //! of the Young's modulus
  //! @see testmaterialspring.h
  PetscScalar getEOmega(void) const { return getE(); }

  //! return Poisson's ratio
  //! @see testmaterialspring.h
  PetscReal getNu(void) const { return 0; }

  //! spring constant for structural interface elements
  //! @see testmaterialspring.h
  PetscScalar getSpringConst(void) const { return 0; }

  //! return the name of the file storing the values
  inline std::string getFilenameCx() const { return m_Filename_Cx; }
  inline std::string getFilenameCy() const { return m_Filename_Cy; }
  inline std::string getFilenameCz() const { return m_Filename_Cz; }
  inline std::string getFilenameCrx() const { return m_Filename_Crx; }
  inline std::string getFilenameCry() const { return m_Filename_Cry; }
  inline std::string getFilenameCrz() const { return m_Filename_Crz; }
  inline std::string getFilenameEtaSpring() const {
    return m_Filename_Eta_Spring;
  }

  //! update the material paramters that are frequency dependent
  void updateMaterial(void);

  /*BEGIN_NO_COVERAGE*/
  //! return the thickness; new: cElement *elemPtr = nullptr
  //! @see testmaterialspring.h
  PetscReal getT(cElement *elemPtr = NULL) const {
    trace(
        "getT() function not defined in isotropic linear elastic spring "
        "material");
    ExitApp();
    return 0.0;
  }
  /*END_NO_COVERAGE*/

  /*BEGIN_NO_COVERAGE*/
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

  bool isnumber(std::string s);
  /*END_NO_COVERAGE*/
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cMaterialSpring &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cMaterialSpring &other) {
  return other.write(os);
}

#endif
