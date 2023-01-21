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

#include "elementfem.h"

cElementFEM::cElementFEM(short NumberOfNodes, short NumberOfDofsPerNode,
                         short NumberOfGaussPoints)
    : cElement(NumberOfNodes),
      m_NumberOfDofsPerNode(NumberOfDofsPerNode),
      m_NumberOfGaussPoints(NumberOfGaussPoints) {
  // empty
}

cElementFEM::cElementFEM(const cElementFEM &other)
    : cElement(other),
      m_NumberOfDofsPerNode(other.getNumberOfDofsPerNode()),
      m_NumberOfGaussPoints(other.getNumberOfGaussPoints()) {
  // --- copy element loads
  for (ItLoadsOnElementMap it = other.getFirstElementLoad();
       it != other.getLastElementLoad(); it++)
    m_ElementLoads.insert(
        std::pair<const short, cElementLoad *>(it->first, it->second));
}

cElementFEM::~cElementFEM() {
  // empty
}

// ---------------------------------------------------------------------------
//   initialize static members
// ---------------------------------------------------------------------------
cMatrix cElementFEM::Jac(3, 3);
PetscReal cElementFEM::detJac = 0.0;
cMatrix cElementFEM::invJac(3, 3);
std::vector<PetscInt> cElementFEM::getGlobalPositions(void) const {
  const int nnod = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();
  std::vector<eKnownDofs> ElementsDofs(getDofs());
  std::vector<PetscInt> rowscols(nnod * ndofs);

  for (int k = 0; k < nnod; k++) {
    for (int dof = 0; dof < (int)ElementsDofs.size(); dof++) {
      rowscols[k * ndofs + dof] = m_Nodes[k]->getGlobalRow(ElementsDofs[dof]);
    }
  }

  return rowscols;
}

void cElementFEM::insertElementLoad(const int &face,
                                    cElementLoad *ptrElementLoad) {
#ifdef PETSC_USE_DEBUG
  if (ptrElementLoad == NULL)
    throw cException("pointer to elementload NULL", __FILE__, __LINE__);
#endif

  m_ElementLoads.insert(
      std::pair<const short, cElementLoad *>(face, ptrElementLoad));
}

void cElementFEM::setupJacobian2D(cArray3d N, int gp) {
  Jac.setValue(0.0);

  // this work for elements in x-y-plane. PlShell4 and PlShell9 include
  // transformation via lambda matrix. Other plate and disc elements don't, so
  // use them solely in x-y-plane.

  for (int k = 0; k < getNumberOfNodes(); k++) {
    // std::cout << "(*m_Nodes[k])[0]: " << (*m_Nodes[k])[0] <<
    // "(*m_Nodes[k])[1]: " << (*m_Nodes[k])[1] << std::endl << std::flush;
    Jac(0, 0) += N(1, k, gp) * (*m_Nodes[k])[0];
    Jac(0, 1) += N(1, k, gp) * (*m_Nodes[k])[1];
    Jac(1, 0) += N(2, k, gp) * (*m_Nodes[k])[0];
    Jac(1, 1) += N(2, k, gp) * (*m_Nodes[k])[1];
  }

  detJac = Jac(0, 0) * Jac(1, 1) - Jac(0, 1) * Jac(1, 0);

  // check orientatiom of elements,
  // if |Jac< < 0.0 something is wrong!
  if (detJac <= 0.0) {
    message("detJac < 0 in element no. %d\n", getId());
    throw cException("check orientation of the elements", __FILE__, __LINE__);
  }
}

void cElementFEM::setupJacobian2D_plate(cArray3d N, int gp,
                                        cArray3d N_map)  // KR
{
  Jac.setValue(0.0);
  int ngp = getNumberOfGaussPoints();

  for (int k = 0; k < getNumberOfNodes(); k++) {
    Jac(0, 0) += N(1, k, gp) * (*m_Nodes[k])[0];
    Jac(0, 1) += N(1, k, gp) * (*m_Nodes[k])[1];
    Jac(1, 0) += N(2, k, gp) * (*m_Nodes[k])[0];
    Jac(1, 1) += N(2, k, gp) * (*m_Nodes[k])[1];
    Jac(2, 0) += 0.;
    Jac(2, 1) += 0.;
  }
  PetscReal g11 =
      Jac(0, 0) * Jac(0, 0) + Jac(1, 0) * Jac(1, 0) + Jac(2, 0) * Jac(2, 0);
  PetscReal g12 =
      Jac(0, 0) * Jac(0, 1) + Jac(1, 0) * Jac(1, 1) + Jac(2, 0) * Jac(2, 1);
  PetscReal g21 = g12;
  PetscReal g22 =
      Jac(0, 1) * Jac(0, 1) + Jac(1, 1) * Jac(1, 1) + Jac(2, 1) * Jac(2, 1);

  detJac = std::sqrt(g11 * g22 - g12 * g21);

  PetscReal inv_det = 1. / (g11 * g22 - g12 * g21);

  // check orientatiom of elements,
  // if |Jac< < 0.0 something is wrong!
  if (detJac <= 0.0) {
    message("detJac < 0 in element no. %d\n", getId());
    throw cException("check orientation of the elements", __FILE__, __LINE__);
  }

  PetscReal g11inv = g22 * inv_det;
  PetscReal g12inv = -g12 * inv_det;
  PetscReal g21inv = -g21 * inv_det;
  PetscReal g22inv = g11 * inv_det;

  PetscReal dxidx_map, dxidy_map, dxidz_map;
  PetscReal detadx_map, detady_map, detadz_map;

  dxidx_map = g11inv * Jac(0, 0) + g12inv * Jac(0, 1);
  dxidy_map = g11inv * Jac(1, 0) + g12inv * Jac(1, 1);
  dxidz_map = g11inv * Jac(2, 0) + g12inv * Jac(2, 1);

  detadx_map = g21inv * Jac(0, 0) + g22inv * Jac(0, 1);
  detady_map = g21inv * Jac(1, 0) + g22inv * Jac(1, 1);
  detadz_map = g21inv * Jac(2, 0) + g22inv * Jac(2, 1);

  for (int k = 0; k < getNumberOfNodes(); k++) {
    N_map(1, k, gp) = N(1, k, gp) * dxidx_map + N(2, k, gp) * detadx_map;
    N_map(2, k, gp) = N(1, k, gp) * dxidy_map + N(2, k, gp) * detady_map;
  }
}

void cElementFEM::setupJacobian3D(cArray3d N, int gp) {
  int nnod = getNumberOfNodes();

  Jac.setValue(0.0);

  // std::cout<<"total number of nodes "<<nnod<<"\n";
  for (int k = 0; k < nnod; k++) {
    Jac(0, 0) += N(1, k, gp) * (*m_Nodes[k])[0];
    Jac(0, 1) += N(1, k, gp) * (*m_Nodes[k])[1];
    Jac(0, 2) += N(1, k, gp) * (*m_Nodes[k])[2];
    Jac(1, 0) += N(2, k, gp) * (*m_Nodes[k])[0];
    Jac(1, 1) += N(2, k, gp) * (*m_Nodes[k])[1];
    Jac(1, 2) += N(2, k, gp) * (*m_Nodes[k])[2];
    Jac(2, 0) += N(3, k, gp) * (*m_Nodes[k])[0];
    Jac(2, 1) += N(3, k, gp) * (*m_Nodes[k])[1];
    Jac(2, 2) += N(3, k, gp) * (*m_Nodes[k])[2];
  }

  detJac =
      Jac(0, 0) * Jac(1, 1) * Jac(2, 2) + Jac(0, 1) * Jac(1, 2) * Jac(2, 0) +
      Jac(0, 2) * Jac(1, 0) * Jac(2, 1) - Jac(0, 2) * Jac(1, 1) * Jac(2, 0) -
      Jac(0, 0) * Jac(1, 2) * Jac(2, 1) - Jac(0, 1) * Jac(1, 0) * Jac(2, 2);

  // check the orientation of the element
  // correct numbering sequence of nodes? If it's wrong then |Jac|<0.0 !
  if (detJac <= 0.0) {
    message("detJac <= 0 in element no. %d\n", getId());
    std::cout << "detJac: " << detJac << "\n";
    // throw cException("check orientation of the elements", __FILE__,
    // __LINE__);
  }
  // else message("detJac ok in element no. %d\n", getId());
  //}
}

void cElementFEM::invertJacobian2D(void) {
  invJac(0, 0) = Jac(1, 1) / detJac;
  invJac(0, 1) = -Jac(0, 1) / detJac;
  invJac(1, 0) = -Jac(1, 0) / detJac;
  invJac(1, 1) = Jac(0, 0) / detJac;
}

void cElementFEM::invertJacobian3D(void) {
  invJac(0, 0) = (Jac(1, 1) * Jac(2, 2) - Jac(1, 2) * Jac(2, 1)) / detJac;
  invJac(1, 0) = -(Jac(1, 0) * Jac(2, 2) - Jac(1, 2) * Jac(2, 0)) / detJac;
  invJac(2, 0) = (Jac(1, 0) * Jac(2, 1) - Jac(1, 1) * Jac(2, 0)) / detJac;
  invJac(0, 1) = -(Jac(0, 1) * Jac(2, 2) - Jac(0, 2) * Jac(2, 1)) / detJac;
  invJac(1, 1) = (Jac(0, 0) * Jac(2, 2) - Jac(0, 2) * Jac(2, 0)) / detJac;
  invJac(2, 1) = -(Jac(0, 0) * Jac(2, 1) - Jac(0, 1) * Jac(2, 0)) / detJac;
  invJac(0, 2) = (Jac(0, 1) * Jac(1, 2) - Jac(0, 2) * Jac(1, 1)) / detJac;
  invJac(1, 2) = -(Jac(0, 0) * Jac(1, 2) - Jac(0, 2) * Jac(1, 0)) / detJac;
  invJac(2, 2) = (Jac(0, 0) * Jac(1, 1) - Jac(0, 1) * Jac(1, 0)) / detJac;
}

std::vector<PetscInt> cElementFEM::getIdsOfFaceNodes(int Face) const {
  std::vector<PetscInt> result(getNumberOfNodesPerFace());
  std::vector<short> indices(getIndicesOfFaceNodes(Face));

  for (int k = 0; k < getNumberOfNodesPerFace(); k++)
    result[k] = getNode(indices[k])->getId();

  return result;
}

void cElementFEM::assembleInternalForces(Vec *fullSolution, cElementMatrix &KM,
                                         cElementVector &LV) {
  // Vector for the local displacements
  cElementVector localSolution(LV.size());
  getLocalSolutionFromFullSolution(fullSolution, localSolution);
  // internal force vector
  cElementVector FV(LV.size());
  infam::mult(KM, localSolution, FV);
  for (int i = 0; i < int(FV.size()); i++) LV[i] -= FV[i];
}

void cElementFEM::getLocalSolutionFromFullSolution(
    Vec *fullSolution, cElementVector &localSolution) {
  if (fullSolution != NULL) {
    // Get Dofs
    int numberOfNodes = this->getNumberOfNodes();
    std::vector<eKnownDofs> dofvector = this->getDofs();
    int k = 0;
    PetscScalar dummy;
    PetscInt index;
    for (int i = 0; i < numberOfNodes; ++i) {
      // std::cout<<"element id "<<this->getId()<<"\n";
      for (int j = 0; j < int(dofvector.size()); ++j) {
        index = this->m_Nodes[i]->getGlobalRow(dofvector.at(j));
        VecGetValues((fullSolution)[0], 1, &index, &dummy);
        localSolution[k] = dummy;
        // std::cout<<"id "<<k<<"\n";
        // std::cout<<"local solution "<<localSolution[k];
        ++k;
      }
    }
  }
}
