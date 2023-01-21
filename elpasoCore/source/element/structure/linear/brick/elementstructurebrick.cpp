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

#include "elementstructurebrick.h"

cElementStructureBrick::cElementStructureBrick(short NumberOfNodes,
                                               short NumberOfGaussPoints)
    : cElementStructureLinear(NumberOfNodes, 3, NumberOfGaussPoints) {
  // empty
}

cElementStructureBrick::cElementStructureBrick(
    const cElementStructureBrick &other)
    : cElementStructureLinear(other) {
  // empty
}

cElementStructureBrick::~cElementStructureBrick() {
  // empty
}

std::vector<eKnownDofs> cElementStructureBrick::getDofs(void) const {
  std::vector<eKnownDofs> res(3);

  res[0] = disp_x1;
  res[1] = disp_x2;
  res[2] = disp_x3;

  return res;
}

void cElementStructureBrick::computeStressesCells(Vec &fullSolution,
                                                  Vec &stresses,
                                                  Vec &stressesSec2) {
  // The stress computation is using a reduced aproach.
  // We need 8 nodes only (corner nodes)

  // const int ngp   = getNumberOfGaussPoints();
  const int nnod = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();
  int nn = 8;  // stress calculation only at 8 corner nodes
  cElementVector localSolution(nnod * ndofs);
  getLocalSolutionFromFullSolution(&fullSolution, localSolution);
  // cMatrix xyz(nnod,3);
  // getLocalCoordinatesOfElementsNodes(xyz);

  cArray3d Nk;  //  shape functions evaluated at nodal points
  Nk.initialize(4, nnod, 1);

  cMatrix Bm(6, nn * ndofs);  // Matrix of Ansatzfunctions

  cElementVector localSolutionHex8(nn * ndofs);  // reduce nodes to 8!

  for (int j = 0; j < nn * ndofs; j++)
    localSolutionHex8[j] =
        localSolution[j];  // coping the local solution form Hex20/27 to Hex8

  cElementVector strain(6);  // local strain
  cElementVector stress(6);  // local stresses

  PetscReal Nx;  // N,x
  PetscReal Ny;  // N,y
  PetscReal Nz;  // N,z
  cPoint ep;

  // -------------------------------------------------------------------------
  //  Elastizitaetsmatrix bestimmen
  // -------------------------------------------------------------------------
  if (m_history != NULL) m_Material->setHistory(m_history);
  cElementMatrix Cm3d(6, 6);
  m_Material->setupC(Cm3d);

  // -- get single node of element. The stress values are accurate at centroid
  // of linear hexahedron element
  ep[0] = 0.;
  ep[1] = 0.;
  ep[2] = 0.;

  for (int f = 0; f < nn; f++) {
    // -- evaluate shape functions at nodes
    Nk(0, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_fun, ep);
    Nk(1, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_xi, ep);
    Nk(2, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_eta, ep);
    Nk(3, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_zeta, ep);
  }

  // -----------------------------------------------------------------------
  //  Aufstellen der Jacobi-Matrix
  // -----------------------------------------------------------------------
  setupJacobian3D(Nk, 0);
  invertJacobian3D();

  // -----------------------------------------------------------------------
  //  Aufstellen der Matrix der Ansatzfkt.
  // -----------------------------------------------------------------------
  for (int k = 0; k < nn; k++) {
    Nx = (invJac(0, 0) * Nk(1, k, 0) + invJac(0, 1) * Nk(2, k, 0) +
          invJac(0, 2) * Nk(3, k, 0));  // dN/dx
    Ny = (invJac(1, 0) * Nk(1, k, 0) + invJac(1, 1) * Nk(2, k, 0) +
          invJac(1, 2) * Nk(3, k, 0));  // dN/dy
    Nz = (invJac(2, 0) * Nk(1, k, 0) + invJac(2, 1) * Nk(2, k, 0) +
          invJac(2, 2) * Nk(3, k, 0));  // dN/dz

    Bm(0, 0 + k * ndofs) = Nx;
    Bm(1, 1 + k * ndofs) = Ny;
    Bm(2, 2 + k * ndofs) = Nz;

    Bm(3, 0 + k * ndofs) = Ny;
    Bm(3, 1 + k * ndofs) = Nx;

    Bm(4, 1 + k * ndofs) = Nz;
    Bm(4, 2 + k * ndofs) = Ny;

    Bm(5, 0 + k * ndofs) = Nz;
    Bm(5, 2 + k * ndofs) = Nx;
  }

  // get local strain
  infam::mult(Bm, localSolutionHex8, strain);
  // get local stress
  infam::mult(Cm3d, strain, stress);

  if (m_history != NULL) {
    // if element has yielding get absolute stress from element history
    if (m_history->getYielding() == true &&
        m_history->hasElasticValues() == true) {
      // std::cout<<"Element "<<getId()<<" getAbsoluteStress(stress)\n";
      m_history->getAbsoluteStress(stress);
      // std::cout<<"stress \n"<<stress<<"\n";
    }
  }
  // add element stresses to global stress vector
  PetscInt pos;

  for (int i = 0; i < 6; ++i) {
    pos = (this->getId0n()) * 9 + i;
    VecSetValue(stresses, pos, stress[i], INSERT_VALUES);
  }
}

void cElementStructureBrick::computeStressesNodes(Vec &fullSolution,
                                                  Vec &stresses) {
  // The stress computation is using a reduced aproach.
  // We need 8 nodes only (corner nodes)

  // const int ngp   = getNumberOfGaussPoints();
  const int nnod = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();
  int nn = 8;  // stress calculation only at 8 corner nodes
  cElementVector localSolution(nnod * ndofs);
  getLocalSolutionFromFullSolution(&fullSolution, localSolution);
  // cMatrix xyz(nnod,3);
  // getLocalCoordinatesOfElementsNodes(xyz);

  cArray3d Nk;  //  shape functions evaluated at nodal points
  Nk.initialize(4, nnod, 1);

  cMatrix Bm(6, nn * ndofs);  // Matrix of Ansatzfunctions

  cElementVector localSolutionHex8(nn * ndofs);  // reduce nodes to 8!

  for (int j = 0; j < nn * ndofs; j++)
    localSolutionHex8[j] =
        localSolution[j];  // coping the local solution form Hex20/27 to Hex8

  cElementVector strain(6);  // local strain
  cElementVector stress(6);  // local stresses

  PetscReal Nx;  // N,x
  PetscReal Ny;  // N,y
  PetscReal Nz;  // N,z
  cPoint ep;

  // -------------------------------------------------------------------------
  //  Elastizitaetsmatrix bestimmen
  // -------------------------------------------------------------------------
  if (m_history != NULL) m_Material->setHistory(m_history);
  cElementMatrix Cm3d(6, 6);
  m_Material->setupC(Cm3d);

  // -- get single node of element. The stress values are accurate at centroid
  // of linear hexahedron element
  ep[0] = 0.;
  ep[1] = 0.;
  ep[2] = 0.;

  for (int f = 0; f < nn; f++) {
    // -- evaluate shape functions at nodes
    Nk(0, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_fun, ep);
    Nk(1, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_xi, ep);
    Nk(2, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_eta, ep);
    Nk(3, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_zeta, ep);
  }

  // -----------------------------------------------------------------------
  //  Aufstellen der Jacobi-Matrix
  // -----------------------------------------------------------------------
  setupJacobian3D(Nk, 0);
  invertJacobian3D();

  // -----------------------------------------------------------------------
  //  Aufstellen der Matrix der Ansatzfkt.
  // -----------------------------------------------------------------------
  for (int k = 0; k < nn; k++) {
    Nx = (invJac(0, 0) * Nk(1, k, 0) + invJac(0, 1) * Nk(2, k, 0) +
          invJac(0, 2) * Nk(3, k, 0));  // dN/dx
    Ny = (invJac(1, 0) * Nk(1, k, 0) + invJac(1, 1) * Nk(2, k, 0) +
          invJac(1, 2) * Nk(3, k, 0));  // dN/dy
    Nz = (invJac(2, 0) * Nk(1, k, 0) + invJac(2, 1) * Nk(2, k, 0) +
          invJac(2, 2) * Nk(3, k, 0));  // dN/dz

    Bm(0, 0 + k * ndofs) = Nx;
    Bm(1, 1 + k * ndofs) = Ny;
    Bm(2, 2 + k * ndofs) = Nz;

    Bm(3, 0 + k * ndofs) = Ny;
    Bm(3, 1 + k * ndofs) = Nx;

    Bm(4, 1 + k * ndofs) = Nz;
    Bm(4, 2 + k * ndofs) = Ny;

    Bm(5, 0 + k * ndofs) = Nz;
    Bm(5, 2 + k * ndofs) = Nx;
  }

  // get local strain
  infam::mult(Bm, localSolutionHex8, strain);
  // get local stress
  infam::mult(Cm3d, strain, stress);

  if (m_history != NULL) {
    // if element has yielding get absolute stress from element history
    if (m_history->getYielding() == true &&
        m_history->hasElasticValues() == true) {
      // std::cout<<"Element "<<getId()<<" getAbsoluteStress(stress)\n";
      m_history->getAbsoluteStress(stress);
      // std::cout<<"stress \n"<<stress<<"\n";
    }
  }
  // add element stresses to global stress vector
  PetscInt pos;

  cElementVector dummy(6);

  for (int k = 0; k < getNumberOfNodes(); k++) {
    // take a copy of stress vector and scale it
    for (int j = 0; j < 6; ++j) dummy[j] = stress[j];
    infam::scale(dummy, double(1. / getNode(k)->getNumEle()));

    for (int i = 0; i < 6; ++i) {
      pos = (getNode(k)->getGlobalSeqId()) * 9 + i;
      VecSetValue(stresses, pos, dummy[i], ADD_VALUES);
    }
  }
}

void cElementStructureBrick::assembleStiffnessMatrix(cElementMatrix &KM,
                                                     Vec *x = NULL,
                                                     Vec *dx = NULL) {
  if (m_history != NULL) {
    m_Material->setHistory(m_history);
  }

  const int ngp = getNumberOfGaussPoints();
  const int nnod = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();
  cMatrix Bm(6, nnod * ndofs);  // Matrix der Ansatzfkt.

  PetscReal Nx;  // N,x
  PetscReal Ny;  // N,y
  PetscReal Nz;  // N,z

  PetscReal wdJ;  // Gewicht * detJ
  cArray3d N(getShapeFunctions());

  // -------------------------------------------------------------------------
  //  Elastizitaetsmatrix bestimmen
  // -------------------------------------------------------------------------
  cElementMatrix Cm3d(6, 6);

  if (getMaterial()->isNonlinearElement() == true) {
    if (m_history == NULL) {
      // std::cout<<"cElementStructureBrick::assembleStiffnessMatrix ->
      // if(getMaterial()->isNonlinearElement()==true) -> init history!\n";
      m_history = new cHistory;
      m_history->setId(getId());
      m_Material->setHistory(m_history);
    }

    // The stress computation is using a reduced aproach.
    // we need 8 nodes only (corner nodes)
    int nn = 8;

    cElementVector localSolution(nnod * ndofs);
    getLocalSolutionFromFullSolution(x, localSolution);

    cArray3d(Nk);  // shape funktions for the reduces approach
    Nk.initialize(4, nnod, 1);

    cMatrix B(6, nn * ndofs);  // Matrix der Ansatzfkt.

    cElementVector localSolutionHex8(nn * ndofs);  // reduced to 8 nodes!

    for (int j = 0; j < nn * ndofs; ++j)
      localSolutionHex8[j] =
          localSolution[j];  // coping the local solution from Hex20/27 to Hex8

    cElementVector strain(6);  // local strain
    cElementVector stress(6);  // local stresses

    // set single node. The stress values are accurate at centroid of linear
    // hexahedron element. This node is not part of Hex8!
    cPoint ep;
    ep[0] = 0.0;
    ep[1] = 0.0;
    ep[2] = 0.0;

    for (int f = 0; f < nn; f++) {
      // evaluate shape funtions at nodes
      Nk(0, f, 0) =
          m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_fun, ep);
      Nk(1, f, 0) =
          m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_xi, ep);
      Nk(2, f, 0) =
          m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_eta, ep);
      Nk(3, f, 0) =
          m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nn, f, N_zeta, ep);
    }

    m_Material->setupC(Cm3d);

    // -----------------------------------------------------------------------
    //  Aufstellen der Jacobi-Matrix
    // -----------------------------------------------------------------------
    setupJacobian3D(Nk, 0);  // has to be 0, due to the reduced approach
    invertJacobian3D();
    // -----------------------------------------------------------------------
    //  Aufstellen der Matrix der Ansatzfkt.
    // -----------------------------------------------------------------------
    for (int k = 0; k < nn; k++) {
      Nx = (invJac(0, 0) * Nk(1, k, 0) + invJac(0, 1) * Nk(2, k, 0) +
            invJac(0, 2) * Nk(3, k, 0));
      Ny = (invJac(1, 0) * Nk(1, k, 0) + invJac(1, 1) * Nk(2, k, 0) +
            invJac(1, 2) * Nk(3, k, 0));
      Nz = (invJac(2, 0) * Nk(1, k, 0) + invJac(2, 1) * Nk(2, k, 0) +
            invJac(2, 2) * Nk(3, k, 0));

      B(0, 0 + k * ndofs) = Nx;
      B(1, 1 + k * ndofs) = Ny;
      B(2, 2 + k * ndofs) = Nz;

      B(3, 0 + k * ndofs) = Ny;
      B(3, 1 + k * ndofs) = Nx;

      B(4, 1 + k * ndofs) = Nz;
      B(4, 2 + k * ndofs) = Ny;

      B(5, 0 + k * ndofs) = Nz;
      B(5, 2 + k * ndofs) = Nx;
    }

    // get local strain
    infam::mult(B, localSolutionHex8, strain);

    m_history->setStrain(strain);

    // get local stress
    infam::mult(Cm3d, strain, stress);

    // if element has yielding get absolute stress from element history
    if (m_history->getYielding() == true &&
        m_history->hasElasticValues() == true) {
      // std::cout<<"Element "<<getId()<<" getAbsoluteStress(stress)\n";
      m_history->getAbsoluteStress(stress);
      // std::cout<<"stress \n"<<stress<<"\n";
    }
    m_history->setStress(stress);
    m_Material->setHistory(m_history);
  }

  m_Material->setupC(Cm3d);

  if (m_Orientation ==
      "user-def")  // Material matrix is given in global system ->
                   // Transformation needed for user-defined coordinate system
  {
    PetscReal lx = m_OrientationVector[0];
    PetscReal ly = m_OrientationVector[1];
    PetscReal lz = m_OrientationVector[2];

    PetscReal rx = m_OrientationVector[3];
    PetscReal ry = m_OrientationVector[4];
    PetscReal rz = m_OrientationVector[5];

    // Calculate the orthogonal third vector for t (the local z direction)
    PetscReal tx = ly * rz - lz * ry;
    PetscReal ty = lz * rx - lx * rz;
    PetscReal tz = lx * ry - ly * rx;

    cMatrix Gb(6, 6);  // Rotation matrix G for 3D
    Gb(0, 0) = lx * lx;
    Gb(0, 1) = ly * ly;
    Gb(0, 2) = lz * lz;
    Gb(0, 3) = lx * ly;
    Gb(0, 4) = lx * lz;
    Gb(0, 5) = ly * lz;
    Gb(1, 0) = rx * rx;
    Gb(1, 1) = ry * ry;
    Gb(1, 2) = rz * rz;
    Gb(1, 3) = rx * ry;
    Gb(1, 4) = rx * rz;
    Gb(1, 5) = ry * rz;
    Gb(2, 0) = tx * tx;
    Gb(2, 1) = ty * ty;
    Gb(2, 2) = tz * tz;
    Gb(2, 3) = lx * ly;
    Gb(2, 4) = lx * lz;
    Gb(2, 5) = ly * lz;
    Gb(3, 0) = 2 * lx * rx;
    Gb(3, 1) = 2 * ly * ry;
    Gb(3, 2) = 2 * lz * rz;
    Gb(3, 3) = lx * ry + ly * rx;
    Gb(3, 4) = lz * rx + lx * rz;
    Gb(3, 5) = ly * rz + lz * ry;
    Gb(4, 0) = 2 * lx * tx;
    Gb(4, 1) = 2 * ly * ty;
    Gb(4, 2) = 2 * lz * tz;
    Gb(4, 3) = tx * ly + ty * lx;
    Gb(4, 4) = tz * lx + tx * lz;
    Gb(4, 5) = ty * lz + tz * ly;
    Gb(5, 0) = 2 * rx * tx;
    Gb(5, 1) = 2 * ry * ty;
    Gb(5, 2) = 2 * rz * tz;
    Gb(5, 3) = rx * ty + ry * tx;
    Gb(5, 4) = rz * tx + rx * tz;
    Gb(5, 5) = ry * tz + rz * ty;

    cElementMatrix Cm3dTemp(6, 6);
    infam::BT_C_B_wdJ(Cm3dTemp, Gb, Cm3d, 1.0);  // Rotation of material matrix
    Cm3d = Cm3dTemp;
  }

  // -------------------------------------------------------------------------
  //  Gauss-Punkt-Schleife
  // -------------------------------------------------------------------------
  for (int n = 0; n < ngp * ngp * ngp; n++) {
    // -----------------------------------------------------------------------
    //  Aufstellen der Jacobi-Matrix
    // -----------------------------------------------------------------------
    setupJacobian3D(N, n);
    invertJacobian3D();

    // -----------------------------------------------------------------------
    //  Aufstellen der Matrix der Ansatzfkt.
    // -----------------------------------------------------------------------
    for (int k = 0; k < nnod; k++) {
      Nx = (invJac(0, 0) * N(1, k, n) + invJac(0, 1) * N(2, k, n) +
            invJac(0, 2) * N(3, k, n));
      Ny = (invJac(1, 0) * N(1, k, n) + invJac(1, 1) * N(2, k, n) +
            invJac(1, 2) * N(3, k, n));
      Nz = (invJac(2, 0) * N(1, k, n) + invJac(2, 1) * N(2, k, n) +
            invJac(2, 2) * N(3, k, n));

      Bm(0, 0 + k * ndofs) = Nx;
      Bm(1, 1 + k * ndofs) = Ny;
      Bm(2, 2 + k * ndofs) = Nz;

      Bm(3, 0 + k * ndofs) = Ny;
      Bm(3, 1 + k * ndofs) = Nx;

      Bm(4, 1 + k * ndofs) = Nz;
      Bm(4, 2 + k * ndofs) = Ny;

      Bm(5, 0 + k * ndofs) = Nz;
      Bm(5, 2 + k * ndofs) = Nx;
    }

    // -----------------------------------------------------------------------
    //  Addieren zur Elementsteifigkeitsmatrix
    // -----------------------------------------------------------------------
    wdJ = m_GaussPoints.getGaussWeight3D(ngp, n) * detJac;

    BT_C_B_wdJ(KM, Bm, Cm3d, wdJ);
  }
}

void cElementStructureBrick::assembleMassMatrix(cElementMatrix &MM) {
  const int ngp = getNumberOfGaussPoints();
  const int nnod = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();

  cMatrix Bm(3, nnod * ndofs);  // Matrix der Ansatzfkt.
  PetscReal wdJ;                // Gewicht fuer akt. Gauss-Pkt * detJ
  cArray3d N(getShapeFunctions());

  // -------------------------------------------------------------------------
  //  benoetigte Materialparameter
  // -------------------------------------------------------------------------
  const PetscReal rho = m_Material->getRho();

  cMatrix Cm(3, 3);
  Cm(0, 0) = rho;
  Cm(1, 1) = rho;
  Cm(2, 2) = rho;

  // -------------------------------------------------------------------------
  //  Gauss-Punkt-Schleife
  // -------------------------------------------------------------------------
  for (int n = 0; n < ngp * ngp * ngp; n++) {
    // -----------------------------------------------------------------------
    //  Aufstellen der Jacobi-Matrix
    // -----------------------------------------------------------------------
    setupJacobian3D(N, n);

    // -----------------------------------------------------------------------
    //  Aufstellen der Matrix der Ansatzfkt.
    // -----------------------------------------------------------------------
    for (int k = 0; k < nnod; k++) {
      Bm(0, 0 + k * ndofs) = N(0, k, n);
      Bm(1, 1 + k * ndofs) = N(0, k, n);
      Bm(2, 2 + k * ndofs) = N(0, k, n);
    }

    // -----------------------------------------------------------------------
    //  Addieren zur Elementmassenmatrix
    // -----------------------------------------------------------------------
    wdJ = m_GaussPoints.getGaussWeight3D(ngp, n) * detJac;
    BT_C_B_wdJ(MM, Bm, Cm, wdJ);
  }
}

void cElementStructureBrick::assembleLoadVector(cElementVector &LV,
                                                cElementMatrix &KM,
                                                Vec *x = NULL, Vec *dx = NULL) {
  // --- leave this function if no loads are applied to this element
  if (m_ElementLoads.size() != 0) {
    std::multimap<short, cElementLoad *>::const_iterator itLoads;
    cArray3d Nface(getShapeFunctionsFace());

    for (itLoads = m_ElementLoads.begin(); itLoads != m_ElementLoads.end();
         itLoads++) {
      // --- check type of load
      cElementLoadStructure *ptrL =
          dynamic_cast<cElementLoadStructure *>(itLoads->second);

      const int nnod_face = getNumberOfNodesPerFace();
      std::vector<short> Inzidenz = getIndicesOfFaceNodes(itLoads->first - 1);
      cElementVector F(3 * nnod_face);
      cMatrix C(nnod_face, nnod_face);

      for (int k = 0; k < nnod_face; k++) {
        F[3 * k] = ptrL->getForceComponent(m_Nodes[Inzidenz[k]]->getId(), 0);
        F[3 * k + 1] =
            ptrL->getForceComponent(m_Nodes[Inzidenz[k]]->getId(), 1);
        F[3 * k + 2] =
            ptrL->getForceComponent(m_Nodes[Inzidenz[k]]->getId(), 2);
      }

      evaluateSurfaceIntegral(itLoads->first - 1, Nface, C);

      for (int z = 0; z < nnod_face; z++) {
        for (int s = 0; s < nnod_face; s++) {
          LV[Inzidenz[z] * 3] += C(z, s) * F[3 * s];
          LV[Inzidenz[z] * 3 + 1] += C(z, s) * F[3 * s + 1];
          LV[Inzidenz[z] * 3 + 2] += C(z, s) * F[3 * s + 2];
        }
      }
    }
  }
  assembleInternalForces(x, KM, LV);
}

void cElementStructureBrick::evaluateSurfaceIntegral(int Face, cArray3d &Nface,
                                                     cMatrix &C) {
  std::vector<short> SurfaceNodes(
      getIndicesOfFaceNodes(Face));  // nodes on one element's face
  std::vector<PetscReal> weights;
  const int nnod_face =
      (int)SurfaceNodes.size();  // number of nodes of one face
  int level;                     // number of Gauss points to use
  int ngps;
  cVector a1(3);  // vectors describing tangential plane
  cVector a2(3);  // ... at current Gauss point

  if (nnod_face == 3) {
    level = 1;
    ngps = level;
    weights.resize(ngps);
    for (int k = 0; k < ngps; k++)
      weights[k] = 0.5 * m_GaussPoints.getGaussWeightTria(level, k);
  } else if (nnod_face == 4) {
    level = 2;
    ngps = level * level;
    weights.resize(ngps);
    for (int k = 0; k < ngps; k++)
      weights[k] = m_GaussPoints.getGaussWeight2D(level, k);
  } else if (nnod_face == 9) {
    level = 3;
    ngps = level * level;
    weights.resize(ngps);
    for (int k = 0; k < ngps; k++)
      weights[k] = m_GaussPoints.getGaussWeight2D(level, k);
  } else {
    // should not get here
    throw cException(
        "Error: do not know what to do for this value of nnod_face !", __FILE__,
        __LINE__);
  }

  for (int n = 0; n < ngps; n++) {
    // --- compute detJac, see e.g.
    //     Boundary Integral Equation Methods for Solids and Fluids
    //     Marc Bonnet
    //     John Wiley & Sons Ltd. (1995)
    //     mech J 111
    a1.setValue(0.0);
    a2.setValue(0.0);

    for (int k = 0; k < nnod_face; k++) {
      for (int d = 0; d < 3; d++) {
        a1[d] += Nface(1, k, n) * (*m_Nodes[SurfaceNodes[k]])[d];
        a2[d] += Nface(2, k, n) * (*m_Nodes[SurfaceNodes[k]])[d];
      }
    }

    detJac = sqrt(a1.abs2() * a2.abs2() - (a1.dot(a2)) * (a1.dot(a2)));

    if (detJac <= 0.0) {
      throw cException(
          "detJac <= 0.0 in cElementFluid::evaluateSurfaceIntegral()", __FILE__,
          __LINE__);
    }

    for (int z = 0; z < nnod_face; z++)
      for (int s = 0; s < nnod_face; s++)
        C(z, s) += Nface(0, z, n) * Nface(0, s, n) * weights[n] * detJac;
  }
}
