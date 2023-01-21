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

#include "elementinterface.h"

cElementInterface::cElementInterface(short NumberOfNodes)
    : cElement(NumberOfNodes), m_MatchingNodes(NumberOfNodes) {
  setRoomId(-1);
  setOrientation(0);

  for (int k = 0; k < NumberOfNodes; k++) m_MatchingNodes[k] = NULL;
}

cElementInterface::cElementInterface(const cElementInterface &other)
    : cElement(other) {
  setRoomId(other.getRoomId());
  setOrientation(other.getOrientation());

  for (int k = 0; k < other.getNumberOfNodes(); k++)
    setMatchingNode(k, other.getMatchingNode(k));
}

cElementInterface::~cElementInterface() {
  for (int k = 0; k < getNumberOfNodes(); k++) m_MatchingNodes[k] = NULL;
}

cNode *cElementInterface::getMatchingNode(int Index) const {
#ifdef PETSC_USE_DEBUG
  return m_MatchingNodes.at(Index);
#else
  return m_MatchingNodes[Index];
#endif
}

void cElementInterface::setMatchingNode(int Index, cNode *ptrNode) {
#ifdef PETSC_USE_DEBUG
  m_MatchingNodes.at(Index) = ptrNode;
#else
  m_MatchingNodes[Index] = ptrNode;
#endif
}

#if 0
void cElementInterface::applyBoundaryConditions(cElementMatrix &C, cElementMatrix &CT)
{
  std::vector<eKnownDofs> DofsStructure( getDofsStructure() );

  const int nnod  = getNumberOfNodes();
  const int ndofs = (int)DofsStructure.size();

  // --- apply nodal boundary conditions
  // --- loop over nodes of element
  PetscInt    pos = 0;  // index of row/column
  PetscScalar val = 0.; // prescribed nodal value

  for (int k=0; k<getNumberOfNodes(); k++)
  {
    // --- loop over boundary conditions applied 
    //     to current structural node
    for (ItNodeBCs itBCs = m_MatchingNodes[k]->getFirstBC(); itBCs != m_MatchingNodes[k]->getLastBC(); itBCs++)
    {
      // --- loop over all degrees of freedom
      for (int dof=0; dof < (int)DofsStructure.size(); dof++)
      {
        // --- check if the dof exists at the node and the dof is fixed by the boundarycondition
        if (itBCs->second->checkIfFixed( DofsStructure[dof] ) == true)
        {
			// NB: wavenumber is currently set to zero. No oblique incidence  as BC
			// on interface elements currently possible!!!!
          val = itBCs->second->getPrescribedValue( DofsStructure[dof], NULL, 0.0 );
          pos = k * ndofs + dof;

          for (int col=0; col<nnod; col++)
          {
            C (pos,col) = 0.;
            CT(col,pos) = 0.;
          }
        }
      }
    }

    // --- loop over boundary conditions applied 
    //     to current fluid node
    for (ItNodeBCs itBCs = m_Nodes[k]->getFirstBC(); itBCs != m_Nodes[k]->getLastBC(); itBCs++)
    {
      // --- check if the dof exists at the node and the dof is fixed by the boundarycondition
      if (itBCs->second->checkIfFixed( fluid ) == true)
      {
		  	// NB: wavenumber is currently set to zero. No oblique incidence  as BC
			// on interface elements currently possible!!!!
        val = itBCs->second->getPrescribedValue( fluid , m_Nodes[k], 0.0);
        pos = k;

        for (int col=0; col<nnod; col++)
        {
          C (col,pos) = 0.;
          CT(pos,col) = 0.;
        }
      }
    }
  }
}
#endif

#if 0
void cElementInterface::addContributionFrequencyDomain(const PetscReal &omega, Mat &K, Vec &F)
{
	// -- rho will be zero for poro coupling elements
  const PetscScalar rho = getRhoF();

  // --- the following two factors will be multiplied to C
  //     and C^T, respectively
  PetscScalar factorC  = -1. * getOrientation();
  PetscScalar factorCT;

  // --- check if symmetric fsi formulation requested
/*  if (cElement::checkSymmetry() == true) { // has to be deleted (meike)
    factorCT = factorC;
  }
  else {
    // -- factors are the same for bounded and unbounded case
    factorC  = -1. * getOrientation();
    factorCT = -1. * getOrientation() * rho * omega * omega;
  }*/
  // for symmetric formulation
  factorCT = factorC;

addContribution(K, K, factorC, factorCT);

}
#endif

#if 0
void cElementInterface::addContributionTimeDomain(Mat &K, Mat &M)
{
  const PetscScalar rho = getRhoF();

  // --- the following two factors will be multiplied to C
  //     and C^T, respectively
  const PetscScalar factorC  = -1. * getOrientation();
  const PetscScalar factorCT =       getOrientation() * rho;

  addContribution(K, M, factorC, factorCT);
}
#endif

std::ostream &cElementInterface::write(std::ostream &os) const {
  os << "id = " << getId() << std::endl;
  os << "fluid nodes" << std::endl;
  for (int k = 0; k < getNumberOfNodes(); k++) os << (*m_Nodes[k]) << std::endl;
  os << "struture nodes" << std::endl;
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << (*m_MatchingNodes[k]) << std::endl;

  os << "ori = " << getOrientation() << std::endl;

  return os;
}
