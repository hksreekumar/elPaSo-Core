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

#ifndef INFAM_ELEMENTINTERFACE_STRUCTURAL_H
#define INFAM_ELEMENTINTERFACE_STRUCTURAL_H

#include "./elementinterface.h"

//! @brief base class for coupling acoustic fluid domains and structural domains
//! @author Dirk Clasen
//! @date 16.08.2006
//! m_Nodes[] holds the nodes of the fluid domain, the nodes of the structural
//! domain that correspond to the fluid domain nodes are stored in a different
//! vector.
//! @todo incorporating non-homogenous boundary conditions
class cElementInterfaceStructural
    : public cElementInterface,
      private tCounter<cElementInterfaceStructural> {
 private:
 protected:
  //! @brief this function performs the computation of C and C^T, multiplies
  //! both by the factors C and C^T and inserts them into the system matrices A
  //! and B by calling insertMatrices(). A and B are both the stiffnessmatrix
  //! for frequency domain analysis. In time domain, A is the stiffness matrix
  //! and B the massmatrix, respectively.
  void addContribution(Mat &A, Mat &B, const PetscScalar &factorC,
                       const PetscScalar &factorCT);

  //! @brief insert the coupling matrices to global system matrices. C is added
  //! to A and C^T is added to B. In frequency domain analysis, A and B are the
  //! same matrix whereas in time domain A is the stiffnessmatrix and B the
  //! massmatrix, respectively.
  void insertMatrices(Mat &A, Mat &B, cElementMatrix &C, cElementMatrix &CT);

  //! @brief insert the coupling matrices to global system matrices. C is added
  //! to A and C^T is added to B. In frequency domain analysis, A and B are the
  //! same matrix whereas in time domain A is the stiffnessmatrix and B the
  //! massmatrix, respectively.
  void insertMatrices_PlaneShellPoro2dUP(Mat &A, Mat &B, cElementMatrix &C,
                                         cElementMatrix &CT);

  //! @brief this function performs the computation of C and C^T, multiplies
  //! both by the factors C and C^T and inserts them into the system matrices A
  //! and B by calling insertMatrices(). A and B are both the stiffnessmatrix
  //! for frequency domain analysis. In time domain, A is the stiffness matrix
  //! and B the massmatrix, respectively.
  void addContribution_PlaneShellPoro2dUP(Mat &A, Mat &B,
                                          const PetscScalar &factorC,
                                          const PetscScalar &factorCT);

  //! @brief Function computes and insert the element load vector into global
  //! load vector. Currently only used for cElementInterfaceStructuralPoro3d
  //! (otherwise coupling does not influence rhs)
  // virtual void addContributionLoadVector(Vec &F) =0;

  //! @brief compute the load vector (rhs) - currently only for coupling fluid
  //! flow
  // virtual void computeLoadVector(cElementVector &LV) = 0;

 public:
  //! @brief Constructor
  cElementInterfaceStructural(short NumberOfNodes);
  //! @brief Copy constructor
  cElementInterfaceStructural(const cElementInterfaceStructural &other);
  //! @brief Destructor
  virtual ~cElementInterfaceStructural();

  //! @brief returns the number of instances of this objecttype
  static size_t howMany(void) {
    return tCounter<cElementInterfaceStructural>::howMany();
  }

  //! @brief apply the nodal boundary conditions to the coupling terms
  void applyBoundaryConditions(cElementMatrix &C, cElementMatrix &CT);

  //! @brief add the contribution of the interfaceelement to the global
  //! matrix and vector. This is the version for the frequency domain analysis.
  //! @param omega angular frequency
  //! @param K global stiffness matrix
  //! @param F global load vector (used when inhomogenous boundary conditions
  //! are applied)
  //! @todo extend this routine to non-homogenous boundary conditions
  void addContributionFrequencyDomain(const PetscReal &omega, Mat &K, Vec &F);

  //! @brief add the contribution of the interfaceelement to the global
  //! matrix and vector. This is the version for the frequency domain analysis
  //! with velocity potential formulation
  //! @param omega angular frequency
  //! @param K global stiffness matrix
  //! @param F global load vector (used when inhomogenous boundary conditions
  //! are applied)
  //! @author Harikrishnan Sreekumar
  void addContributionFrequencyDomainVelocityPotential(const PetscReal &omega,
                                                       Mat &A, Mat &B, Vec &F);

  //! @brief add the contribution of the interfaceelement to the global
  //! matrix and vector. This is the version for the time domain/eigenvalue
  //! analysis.
  //! @param K global stiffness matrix
  //! @param M global mass matrix
  //! @todo extend this routine to non-homogenous boundary conditions
  void addContributionTimeDomain(Mat &K, Mat &M);

  //! @brief write this object to a stream
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief write this object in XML format to a stream
  //! @param os outputstream
  //! @return der modified outputstream
  // virtual std::ostream& writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &write(std::ostream &os,
                           const cElementInterfaceStructural &other) {
  return other.write(os);
}

#endif
