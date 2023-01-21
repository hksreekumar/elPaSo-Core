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

#ifndef INFAM_ANALYSIS_FREQUENCY_H
#define INFAM_ANALYSIS_FREQUENCY_H

#include "../../element/structure/linear/plate/elementstructurekirchhoff.h"
#include "../analysis.h"

//! @brief computes frequency response - basic routine for elpasoCore
//! @author Dirk Clasen, Harikrishnan Sreekumar
//! @date 15.08.2005
class cAnalysisFrequencyBasic : public cAnalysis {
 private:
  //! first frequency
  PetscReal m_Start;
  //! number of frequencies to compute
  PetscInt m_Steps;
  //! increment (delta f)
  PetscReal m_Increment;
  //! all frequencies to compute
  std::set<PetscReal> m_Frequencies;
  //! iterator for looping frequencies
  std::set<PetscReal>::iterator m_ItFrequencies;

  //! @brief Discretizes the frequency domain
  //! @param Filename The filename to read the freq data
  void generateFrequencies(const std::string &Filename);

 protected:
  //! @brief function to assemble global tensors (Stiffnessmatrix, Massmatrix,
  //! Loadvector)
  //! @param myMesh used mesh
  //! @param nmodes used by similar sv4
  //! @param omega angular freq.
  void assembleGlobalTensors(cMesh &MyMesh, int nmodes = 0,
                             const PetscReal omega = 0);

  //! @brief allocates memory and initializes the solvers.
  //! @param myMesh used mesh
  void initializePETScObjects(cMesh &myMesh);

  //! @brief frees memory and destroys the solver objects.
  void deletePETScObjects(void);

 public:
  //! @brief Constructor
  cAnalysisFrequencyBasic();
  //! @brief Destructor
  ~cAnalysisFrequencyBasic();

  //! @brief read first frequency
  //! @return start frequency
  inline PetscReal getStartFrequency(void) const { return m_Start; }

  //! @brief returns the number of frequencies to compute
  //! @return frequency steps
  inline PetscInt getNumberOfSteps(void) const { return m_Steps; }

  //! @brief returns the increment \f$ (\Delta f) \f$
  //! @return frequency step size
  inline PetscReal getIncrement(void) const { return m_Increment; }

  //! @brief performs a full computation run
  //! @param MyMesh mesh used for computation
  //! @param freemem tells us whether the Petsc Objects are going to be deletet
  //! at the and of FullRun. The default is TRUE
  void FullRun(cMesh &MyMesh, PetscBool freemem = PETSC_TRUE);

  //! @brief perform full run with similar coupling (Carolin Birk)
  //! @param MyMesh mesh to compute on
  //! @param freemem tells us whether the Petsc Objects are going to be deletet
  //! at the and of FullRun. The default is TRUE
  void FullRunSimilarCB(cMesh &MyMesh, PetscBool freemem = PETSC_TRUE);

  //! @brief writes this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief writes this object in XML format to a stream
  //! @brief This function is used after searching for interface elements when a
  //! new inputfile is generated that contains these interfaceelements
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &writeXml(std::ostream &os) const;

  //! @brief reads this object of a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);
};

//! @brief overloaded outputoperator
inline std::istream &operator>>(std::istream &is,
                                cAnalysisFrequencyBasic &other) {
  return other.read(is);
}

//! @brief overloaded inputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cAnalysisFrequencyBasic &other) {
  return other.write(os);
}

#endif
