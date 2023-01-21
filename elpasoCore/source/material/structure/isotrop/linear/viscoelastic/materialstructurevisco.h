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

#ifndef INFAM_MATERIALSTRVISCO_H
#define INFAM_MATERIALSTRVISCO_H

#include <string>

#include "../materialstructureisolin.h"

//! @brief viscoelastic material
//! @author Dirk Clasen, Christopher Blech
//! @date 04.03.2004, revised 07.10.2017
class cMaterialStructureVisco : public cMaterialStructureIsoLin,
                                private tCounter<cMaterialStructureVisco> {
 private:
  std::string m_Filename_Eta;  ///< name of the inputfile for eta
  PetscInt m_lines_in_eta_file;
  std::vector<PetscReal> m_Omega_Eta;  ///< frequencies for eta
  std::vector<PetscReal> m_Eta_var;    ///< Eta, frequency-dependent

  std::string m_Filename_E;  ///< name of the inputfile for E
  PetscInt m_lines_in_E_file;
  std::vector<PetscReal> m_Omega_E;  ///< frequencies for E
  std::vector<PetscReal> m_E_var;    ///< E, frequency-dependent

  std::string m_Filename_T;  ///< name of the inputfile for T
  PetscInt m_lines_in_T_file;
  std::vector<PetscReal>
      m_x_Tdata;  ///< x data from coordinate based thickness file
  std::vector<PetscReal>
      m_y_Tdata;  ///< y data from coordinate based thickness file
  std::vector<PetscReal>
      m_T_var;  ///< thickness data from coordinate based thickness file

  mutable std::vector<std::pair<cElement *, PetscReal> >
      m_calculatedElementTData;  ///< Stores elements with their previously
                                 ///< calculated thickness data
  mutable PetscInt
      m_calculatedElementTDataLength;  ///< Stores length of
                                       ///< m_calculatedElementTData vector

  PetscScalar m_E;    ///< Young's modulus
  PetscReal m_Eta;    ///< damping parameter
  PetscReal m_Nu;     ///< Poisson's ratio
  short m_ViscoType;  ///< flag that defines the type of damping

  //! read file storing the loss factors
  void readEtaFromFile(void);

  //! read file storing the Young's moduli
  void readEFromFile(void);

  //! read file storing the thickness in dependency on the spatial coordinate
  void readTFromFile(void);

 protected:
  PetscScalar E_vis;  ///< viscoelastic Young's modulus (computed for current
                      ///< frequency omega)
  PetscReal m_Ks;     ///< shear correction factor

  std::complex<PetscReal> m_ImagOne;  ///< complex number i=\sqrt{1}

  PetscScalar getEVisco(void) const { return E_vis; }

 public:
  cMaterialStructureVisco();
  cMaterialStructureVisco(const cMaterialStructureVisco &other);
  ~cMaterialStructureVisco();

  //! return the number of instances of this object
  static size_t howMany(void) {
    return tCounter<cMaterialStructureVisco>::howMany();
  }

  //! return frequency dependent Young's modulus
  //! @see testmaterialstructurevisco.h
  PetscScalar getEOmega(void) const { return E_vis; }

  //! @see testmaterialstructurevisco.h
  inline PetscScalar getE(void) const { return m_E; }

  //! @see testmaterialstructurevisco.h
  inline PetscReal getEta() const { return m_Eta; }

  //! @see testmaterialstructurevisco.h
  inline PetscReal getNu() const { return m_Nu; }

  //! @see testmaterialstructurevisco.h
  inline short getType() const { return m_ViscoType; }

  //! @see testmaterialstructurevisco.h
  inline void setE(const PetscScalar &value) { m_E = value; }

  //! @see testmaterialstructurevisco.h
  inline void setEta(const PetscReal &value) { m_Eta = value; }

  //! @see testmaterialstructurevisco.h
  inline void setNu(const PetscReal &value) { m_Nu = value; }

  //! @see testmaterialstructurevisco.h
  inline void setType(const short &value) { m_ViscoType = value; }

  //! return the name of the file storing the values
  //! @see testmaterialstructurevisco.h
  inline std::string getFilenameEta() const { return m_Filename_Eta; }

  //! return the name of the file storing the values
  //! @see testmaterialstructurevisco.h
  inline std::string getFilenameE() const { return m_Filename_E; }

  /*BEGIN_NO_COVERAGE*/
  //! return the name of the file storing the values
  inline std::string getFilenameT() const { return m_Filename_T; }
  /*END_NO_COVERAGE*/

  //! update the material paramters that are frequency dependent
  void updateMaterial(void);

  //! spring constant for structural interface elements
  //! @see testmaterialstructurevisco.h
  PetscScalar getSpringConst(void) const { return getEVisco(); }

  //! get thickness of thickness file according to corresponding element
  //! default value NULL
  //! @see testmaterialstructurevisco.h
  PetscReal getT(cElement *elemPtr = NULL) const;

  //! method to calculate the mean coordinate of a single component (x=0, y=1 or
  //! z=2) of all nodes of an element
  //! @see testmaterialstructurevisco.h
  PetscReal calculateMeanCoordinate(cElement *elementPtr,
                                    int coord_index) const;

  //! method to calculate the distance of two points (x1,y1), (x2,y2) in
  //! xy-plane
  //! @see testmaterialstructurevisco.h
  PetscReal calculateDistance(PetscReal x1, PetscReal y1, PetscReal x2,
                              PetscReal y2) const;

  //! shear part of the elasticity matrix of Mindlin plate
  //! @see testmaterialstructurevisco.h
  void setupCs(cElementMatrix &Cs, cElement *elemPtr = NULL) const;

  //! bending part of the elasticity matrix of Mindlin plate
  //! @see testmaterialstructurevisco.h
  void setupCb(cElementMatrix &Cb, cElement *elemPtr = NULL) const;

  /*BEGIN_NO_COVERAGE*/
  //! prestress part of the elasticity matrix of Mindlin plate
  void setupCpre(cElementMatrix &Cpre, cElement *elemPtr = NULL) const {}
  /*END_NO_COVERAGE*/

  //! elasticity matrix of membrane elements
  //! @see testmaterialstructurevisco.h
  void setupCm(cElementMatrix &Cm) const;

  /*BEGIN_NO_COVERAGE*/
  //! elasticity matrix of 3d solid elements
  void setupC(cElementMatrix &C) const;

  //! elasticity matrix for plane stress problem
  void setupCps(cElementMatrix &Cps) const {}
  /*END_NO_COVERAGE*/

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os) const;

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is,
                                cMaterialStructureVisco &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cMaterialStructureVisco &other) {
  return other.write(os);
}

#endif
