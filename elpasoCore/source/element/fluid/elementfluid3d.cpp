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

#include "elementfluid3d.h"

cElementFluid3d::cElementFluid3d(short NumberOfNodes, short NumberOfGaussPoints)
    : cElementFluid(NumberOfNodes, NumberOfGaussPoints) {
  // empty
}

cElementFluid3d::cElementFluid3d(const cElementFluid3d &other)
    : cElementFluid(other) {
  // empty
}

cElementFluid3d::~cElementFluid3d() {
  // empty
}

cVector cElementFluid3d::getNormalOnFace(const short &face) const {
  cVector normal(3), a(3), b(3);
  std::vector<short> FaceIncidence(getIndicesOfFaceNodes(face));

  for (int d = 0; d < 3; d++) {
    a[d] = (*m_Nodes[FaceIncidence[1]])[d] - (*m_Nodes[FaceIncidence[0]])[d];
    b[d] = (*m_Nodes[FaceIncidence[2]])[d] - (*m_Nodes[FaceIncidence[0]])[d];
  }

  normal[0] = a[1] * b[2] - a[2] * b[1];
  normal[1] = -a[0] * b[2] + a[2] * b[0];
  normal[2] = a[0] * b[1] - a[1] * b[0];

  infam::scale(normal, 1. / normal.abs());

  return normal;
}

void cElementFluid3d::assembleMassMatrix(cElementMatrix &MM) {
  const int nnod = getNumberOfNodes();
  const int ngp = getNumberOfGaussPoints();

  const PetscScalar cf2 =
      1.0 / (m_Material->getCfOmega() * m_Material->getCfOmega());

  PetscReal wdJ;  // cur. Gauss weight multyplied by detJac
  cArray3d N(getShapeFunctions());
  PetscReal rho_c, lambda;
  cPoint cgp;  // current gp

  for (int n = 0; n < ngp * ngp * ngp; n++) {
    setupJacobian3D(N, n);
    cgp[0] = 0.0;
    cgp[1] = 0.0;
    cgp[2] = 0.0;

    wdJ = m_GaussPoints.getGaussWeight3D(ngp, n) * detJac;

    // insert components if cloaking material
    /*  cMaterialFluidCloaking *ptrC = dynamic_cast<cMaterialFluidCloaking
      *>(m_Material); if (ptrC != NULL)
          {
                  // compute real coordinates of GP (current GP is defined by n)
                  for (int k=0; k<nnod; k++)
                  {
                          cgp[0]+=N(0,k,n) * (*getNode(k))[0]   ;
                          cgp[1]+=N(0,k,n) * (*getNode(k))[1];
                          cgp[2]+=N(0,k,n) * (*getNode(k))[2];
                  }
                  // get rho_c for a point and get lambda
                  lambda = getMaterial()-> getRhoLambda (cgp, rho_c);

                  //contribution of anisotropic cloaking material (divide by
      lambda contribution)
                  // Cummer et al. 2008
                  wdJ *=  1 / lambda;
          }*/

    for (int z = 0; z < nnod; z++)
      for (int s = 0; s < nnod; s++)
        MM(z, s) += cf2 * N(0, z, n) * N(0, s, n) * wdJ;
  }
}

void cElementFluid3d::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL,
                                              Vec *dx = NULL) {
  const int nnod = getNumberOfNodes();
  PetscReal wdJ;            // weight * detJac
  cMatrix H(3, nnod);       // [dN/dxi, dN/deta, dN/dzeta]^T at each node
  cMatrix JH(3, nnod);      // = invJ * H
  cMatrix tJH(nnod, 3);     // = JH^T
  cMatrix JtJ(nnod, nnod);  // = tJH * JH
  cArray3d N(getShapeFunctions());
  const int ngp(getNumberOfGaussPoints());
  cPoint cgp;

  for (int n = 0; n < ngp * ngp * ngp; n++) {
    setupJacobian3D(N, n);
    invertJacobian3D();

    // ------------------------------------------------------------------------
    //  matrix H used to compute the gradient:  J^{-1} * H = grad(N)
    // ------------------------------------------------------------------------
    for (int k = 0; k < nnod; k++) {
      H(0, k) = N(1, k, n);
      H(1, k) = N(2, k, n);
      H(2, k) = N(3, k, n);
    }

    wdJ = m_GaussPoints.getGaussWeight3D(ngp, n) * detJac;

    JH.setValue(0.0);
    infam::mult(invJac, H, JH);
    infam::transpose(JH, tJH);

    // insert components if cloaking material
    /*cMaterialFluidCloaking *ptrC = dynamic_cast<cMaterialFluidCloaking
    *>(m_Material); if (ptrC != NULL)
        {
                cPoint cgp;					// current gp
        cMatrix rho_comp(3, 3);		// matrix for components of rho
                cMatrix help(nnod,3);		// auxiliary matrix

                // compute real coordinates of GP (current GP is defined by n)
                for (int k=0; k<nnod; k++)
                {
                        cgp[0]+=N(0,k,n) * (*getNode(k))[0];
                        cgp[1]+=N(0,k,n) * (*getNode(k))[1];
                        cgp[2]+=N(0,k,n) * (*getNode(k))[2];
                }
                // get rho_x, rho_y and rho_z components
                getMaterial()-> computeRhoValues (cgp, rho_comp);
                //rho_comp.setValue(0.0);
                // set rho components in matrix
                //rho_comp(0,0) = rx; rho_comp(1,1) = ry; rho_comp(2,2) = rz;

                infam::mult(tJH, rho_comp, help);
                tJH = help;
        }*/

    cMatrix result(KM.rows(), KM.cols());
    infam::mult(tJH, JH, result);
    infam::scale(result, wdJ);
    infam::add(result, KM);
  }

  // additional factor, if fluid is an equivalent fluid
  /*cMaterialFluidEquiv *ptrE = dynamic_cast<cMaterialFluidEquiv *>(m_Material);
  if (ptrE != NULL) {
        infam::scale(KM, ptrE->getCouplingFactor());
  }*/
}

void cElementFluid3d::evaluateSurfaceIntegral(int Face, cArray3d &Nface,
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

void cElementFluid3d::assembleLoadVector(cElementVector &LV, cElementMatrix &KM,
                                         Vec *x = NULL, Vec *dx = NULL) {
  // empty
}

void cElementFluid3d::assembleDynamicLoadVector(const PetscReal &omega,
                                                cElementVector &LV,
                                                cElementMatrix &EM) {
  // --- leave this function if no loads are applied to this element
  if (m_ElementLoads.size() == 0) return;

  std::multimap<short, cElementLoad *>::const_iterator itLoads;
  cArray3d Nface(getShapeFunctionsFace());

  for (itLoads = m_ElementLoads.begin(); itLoads != m_ElementLoads.end();
       itLoads++) {
    // -----------------------------------------------------------------------
    //   check if load is a normal velocity
    // -----------------------------------------------------------------------
    /*cElementLoadVn *ptrVn = dynamic_cast<cElementLoadVn *>(itLoads->second);
    if (ptrVn != NULL)
    {
      const int          nnod_face = getNumberOfNodesPerFace();
      std::vector<short> Inzidenz  = getIndicesOfFaceNodes(itLoads->first - 1);
      cElementVector     Vn( nnod_face );
      cMatrix            C( nnod_face, nnod_face );

#ifdef PETSC_USE_COMPLEX
      if (ptrVn->getType() == vn)
                    //Vn.setValue( 1. * std::complex<PetscReal>(.0, 1.0) *
m_Material->getRhoOmega() * omega * ptrVn->getValue() ); // by Dirk Vn.setValue(
-1. * std::complex<PetscReal>(.0, 1.0) * m_Material->getRhoOmega() * omega *
ptrVn->getValue() ); // by Meike else if (ptrVn->getType() == flux) Vn.setValue(
ptrVn->getValue()); else
      {
        throw cException("no type specified for cElementLoadVn", __FILE__,
__LINE__);
      }
#else
// vn without complex numbers: assuming time-domain computation
         // if (m_AnalysisType == Time)
        //	{
                        if (ptrVn->getType() == vn || ptrVn->getType() == flux)
                                Vn.setValue( ptrVn->getValue() );
                        else
                        {
                                throw cException("no type specified for
cElementLoadVn", __FILE__, __LINE__);
                        }
        //	}
        //  else
        //	{
        //	throw cException("UNABLE TO EVALUATE ELEMENT LOAD ON FLUID -
NEED COMPLEX NUMBERS", __FILE__, __LINE__);
        //	}
#endif
      evaluateSurfaceIntegral(itLoads->first - 1, Nface, C);

          for (int z=0; z<nnod_face; z++)
      for (int s=0; s<nnod_face; s++)
        LV[Inzidenz[z]] += C(z,s) * Vn[s];
    }
*/
  }

  // --- for symmetric formulation
  const PetscScalar factor = 1. / (m_Material->getRhoOmega() * omega * omega);
  infam::scale(LV, factor);

  // additional factor, if fluid is an equivalent fluid
  /*cMaterialFluidEquiv *ptrE = dynamic_cast<cMaterialFluidEquiv *>(m_Material);
  if (ptrE != NULL) {
        infam::scale(LV, ptrE->getCouplingFactor());
  }*/
}

void cElementFluid3d::evaluateImpedance(const PetscReal &omega,
                                        cElementMatrix &KM) {
  // --- leave this function if no elementloads are applied
  if (m_ElementLoads.size() == 0) return;

  std::multimap<short, cElementLoad *>::const_iterator itLoads;
  cArray3d Nface(getShapeFunctionsFace());

  for (itLoads = m_ElementLoads.begin(); itLoads != m_ElementLoads.end();
       itLoads++) {
    // -----------------------------------------------------------------------
    //   check if elementload is an impedance
    // -----------------------------------------------------------------------
    /*cElementLoadImpedance *ptrZ = dynamic_cast<cElementLoadImpedance
*>(itLoads->second);

    if (ptrZ != NULL)
    {
      // --- get nodes that describe elements' surface
      //     we need the index of the nodes within incidence table
      std::vector<short>       Inzidenz  = getIndicesOfFaceNodes(itLoads->first
- 1); const int                nnod_face = getNumberOfNodesPerFace();
      std::vector<PetscScalar> Z( nnod_face, 0.);
      std::vector<PetscScalar> HD( nnod_face, 0. );
      cMatrix                  C( nnod_face, nnod_face );

#ifdef PETSC_USE_COMPLEX
      if (omega != 0.0)
      {
        // --- compute
        //                 i * rho * omega
        //     factor = ---------------------
        //                        Z
        //
                // if Z = normalized Wallimpedance we have to modify the factor:
                //
        //                 i * rho * omega
        //     factor = ---------------------
        //                   Z * rho * c
        //

        const PetscScalar factor = omega * m_Material->getRhoOmega() *
std::complex<double>(0.,1.) /  ptrZ->getImpedance(omega);

        for (int k=0; k<nnod_face; k++)
          Z[k] = factor;
      }
#endif

      // --- evaluate surfaceintegral
      evaluateSurfaceIntegral(itLoads->first - 1, Nface, C);

      for (int z=0; z<nnod_face; z++)
        for (int s=0; s<nnod_face; s++)
          HD[z] += C(z,s) * Z[s];


      for (int z=0; z<nnod_face; z++)
        KM(Inzidenz[z],Inzidenz[z]) += HD[z];
    }*/
  }

  // additional factor, if fluid is an equivalent fluid
  /*cMaterialFluidEquiv *ptrE = dynamic_cast<cMaterialFluidEquiv *>(m_Material);
  if (ptrE != NULL) {
        infam::scale(KM, ptrE->getCouplingFactor());
  }*/
}

void cElementFluid3d::computeStressesNodes(Vec &FullSolution, Vec &Stresses) {
  // --- here we will compute grad(p) and average it.
  const int nnod = getNumberOfNodes();
  cElementVector gradP(3 * nnod);
  cElementVector valP(nnod);
  cElementVector ones(3 * nnod);
  PetscInt index;

  // --- extract fluid pressure from solution vector
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);

    VecGetValues(FullSolution, 1, &index, &(valP[k]));
  }

  computeGradP(gradP, valP);

  // --- add to global vectors
  PetscInt *indices = new PetscInt[3 * nnod];

  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalSeqId();
    indices[3 * k] = cstNumberOfStressDofs * index + gradpx;
    indices[3 * k + 1] = cstNumberOfStressDofs * index + gradpy;
    indices[3 * k + 2] = cstNumberOfStressDofs * index + gradpz;

    ones[3 * k] = 1.;
    ones[3 * k + 1] = 1.;
    ones[3 * k + 2] = 1.;
  }

  VecSetValues(Stresses, 3 * nnod, indices, &(gradP[0]), ADD_VALUES);

  delete[] indices;
}

void cElementFluid3d::computeGradP(cElementVector &gradP,
                                   const cElementVector &valP) {
  cArray3d Nk;  //  shape functions evaluated at nodal points
  const int nnod = getNumberOfNodes();
  cPoint ep;

  cMatrix xyz(nnod, 3);

  // -- local coordinates of element's nodes
  getLocalCoordinatesOfElementsNodes(xyz);

  Nk.initialize(4, nnod, 1);

  // -----------------------------------------------------------------------
  //  Aufstellen der Matrix der Ansatzfkt.
  // -----------------------------------------------------------------------
  for (int k = 0; k < nnod; k++) {
    cElementVector gradPPoint(3);
    // -- get single node of element
    for (int i = 0; i < 3; i++) {
      ep[i] = xyz(k, i);
    }

    computeGradPPoint(ep, Nk, gradPPoint, valP);

    gradP[3 * k] = gradPPoint[0];
    gradP[3 * k + 1] = gradPPoint[1];
    gradP[3 * k + 2] = gradPPoint[2];
  }
}

void cElementFluid3d::computeGradPPoint(cPoint &ep, cArray3d &Nk,
                                        cElementVector &gradPPoint,
                                        const cElementVector &valP) {
  const int nnod = getNumberOfNodes();
  cMatrix JH(3, nnod);  // [N,x N,y N,z]^T for each node
  PetscScalar valtemp;

  std::vector<PetscReal> gradPPoint_list_Xr(nnod);
  std::vector<PetscReal> gradPPoint_list_Xi(nnod);
  std::vector<PetscReal> gradPPoint_list_Yr(nnod);
  std::vector<PetscReal> gradPPoint_list_Yi(nnod);
  std::vector<PetscReal> gradPPoint_list_Zr(nnod);
  std::vector<PetscReal> gradPPoint_list_Zi(nnod);
  std::vector<PetscReal> sumgradPPoint_list(6);

  gradPPoint.setValue(0.0);

  //	for (int i=0; i<6; i++)
  //		gradPPoint_list[i].resize(nnod);

  // -----------------------------------------------------------------------
  //  Aufstellen der Jacobi-Matrix
  // -----------------------------------------------------------------------
  for (int f = 0; f < nnod; f++) {
    // -- evaluate shape functions at nodes
    Nk(0, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, f, N_fun, ep);
    Nk(1, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, f, N_xi, ep);
    Nk(2, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, f, N_eta, ep);
    Nk(3, f, 0) =
        m_ShapeFunctions.evaluateShapeFunction(Hexahedron, nnod, f, N_zeta, ep);
  }

  // -- compute Jacobian
  setupJacobian3D(Nk, 0);
  invertJacobian3D();

  for (int f = 0; f < nnod; f++) {
    JH(0, f) = (invJac(0, 0) * Nk(1, f, 0) + invJac(0, 1) * Nk(2, f, 0) +
                invJac(0, 2) * Nk(3, f, 0));  // dN/dx
    JH(1, f) = (invJac(1, 0) * Nk(1, f, 0) + invJac(1, 1) * Nk(2, f, 0) +
                invJac(1, 2) * Nk(3, f, 0));  // dN/dy
    JH(2, f) = (invJac(2, 0) * Nk(1, f, 0) + invJac(2, 1) * Nk(2, f, 0) +
                invJac(2, 2) * Nk(3, f, 0));  // dN/dz

    // store the imag and real component of the grad value in a list
    // later the list will be sortet and summed in order to minimize numerical
    // errors

#ifdef PETSC_USE_COMPLEX
    valtemp = JH(0, f) * valP[f];
    gradPPoint_list_Xr[f] = valtemp.real();
    gradPPoint_list_Xi[f] = valtemp.imag();
    valtemp = JH(1, f) * valP[f];
    gradPPoint_list_Yr[f] = valtemp.real();
    gradPPoint_list_Yi[f] = valtemp.imag();
    valtemp = JH(2, f) * valP[f];
    gradPPoint_list_Zr[f] = valtemp.real();
    gradPPoint_list_Zi[f] = valtemp.imag();
#else
    valtemp = JH(0, f) * valP[f];
    gradPPoint_list_Xr[f] = valtemp;
    gradPPoint_list_Xi[f] = 0;
    valtemp = JH(1, f) * valP[f];
    gradPPoint_list_Yr[f] = valtemp;
    gradPPoint_list_Yi[f] = 0;
    valtemp = JH(2, f) * valP[f];
    gradPPoint_list_Zr[f] = valtemp;
    gradPPoint_list_Zi[f] = 0;
#endif
  }

  std::sort(gradPPoint_list_Xr.begin(), gradPPoint_list_Xr.end(), abs_order);
  std::sort(gradPPoint_list_Xi.begin(), gradPPoint_list_Xi.end(), abs_order);
  std::sort(gradPPoint_list_Yr.begin(), gradPPoint_list_Yr.end(), abs_order);
  std::sort(gradPPoint_list_Yi.begin(), gradPPoint_list_Yi.end(), abs_order);
  std::sort(gradPPoint_list_Zr.begin(), gradPPoint_list_Zr.end(), abs_order);
  std::sort(gradPPoint_list_Zi.begin(), gradPPoint_list_Zi.end(), abs_order);

  for (int f = gradPPoint_list_Xr.size() - 1; f >= 0; f--) {
    sumgradPPoint_list[0] += gradPPoint_list_Xr[f];
  }
  for (int f = gradPPoint_list_Xi.size() - 1; f >= 0; f--) {
    sumgradPPoint_list[1] += gradPPoint_list_Xi[f];
  }
  for (int f = gradPPoint_list_Yr.size() - 1; f >= 0; f--) {
    sumgradPPoint_list[2] += gradPPoint_list_Yr[f];
  }
  for (int f = gradPPoint_list_Yi.size() - 1; f >= 0; f--) {
    sumgradPPoint_list[3] += gradPPoint_list_Yi[f];
  }
  for (int f = gradPPoint_list_Zr.size() - 1; f >= 0; f--) {
    sumgradPPoint_list[4] += gradPPoint_list_Zr[f];
  }
  for (int f = gradPPoint_list_Zi.size() - 1; f >= 0; f--) {
    sumgradPPoint_list[5] += gradPPoint_list_Zi[f];
  }

#ifdef PETSC_USE_COMPLEX
  for (int i = 0; i < 3; i++) {
    gradPPoint[i] =
        sumgradPPoint_list[2 * i] +
        std::complex<PetscReal>(.0, 1.0) * sumgradPPoint_list[2 * i + 1];
  }
#else
  for (int i = 0; i < 3; i++) {
    gradPPoint[i] = sumgradPPoint_list[2 * i];
  }
#endif
}

void cElementFluid3d::getGaussPointsFace(int Face, std::vector<cPoint> &gp,
                                         int &ngps,
                                         std::vector<PetscReal> &weights,
                                         const int nnod_face) {
  int level;  // number of Gauss points to use

  if (nnod_face == 3)  //"Error: Not yet implemented for tetra-elements!"
  {
    level = 1;
    ngps = level;
    weights.resize(ngps);
    gp.resize(ngps);
    for (int n = 0; n < ngps; n++) {
      weights[n] = 0.5 * m_GaussPoints.getGaussWeightTria(level, n);
      gp[n] = m_GaussPoints.getGaussPointTria(level, n);
    }
  } else if (nnod_face == 4) {
    level = 2;
    ngps = level * level;
    weights.resize(ngps);
    gp.resize(ngps);
    for (int n = 0; n < ngps; n++) {
      weights[n] = m_GaussPoints.getGaussWeight2D(level, n);
      gp[n] = m_GaussPoints.getGaussPoint2D(level, n);
    }
  } else if (nnod_face == 9) {
    level = 3;
    ngps = level * level;
    weights.resize(ngps);
    gp.resize(ngps);
    for (int n = 0; n < ngps; n++) {
      weights[n] = m_GaussPoints.getGaussWeight2D(level, n);
      gp[n] = m_GaussPoints.getGaussPoint2D(level, n);
    }
  } else {
    // should not get here
    throw cException(
        "Error: do not know what to do for this value of nnod_face !", __FILE__,
        __LINE__);
  }

  if (nnod_face == 4 || nnod_face == 9) {
    for (int n = 0; n < ngps; n++) {
      PetscReal xi_new = 0;
      PetscReal eta_new = 0;
      PetscReal zeta_new = 0;

      switch (Face) {
        case 0:
          // face normal in neative eta direction
          xi_new = (gp[n])[0];
          eta_new = -1.0;
          zeta_new = (gp[n])[1];
          break;
        case 1:
          // face normal in positive eta direction
          xi_new = (gp[n])[0];
          eta_new = 1.0;
          zeta_new = (gp[n])[1];
          break;
        case 2:
          // face normal in neative zeta direction
          xi_new = (gp[n])[0];
          eta_new = (gp[n])[1];
          zeta_new = -1.0;
          break;
        case 3:
          // face normal in positive xi direction
          xi_new = 1.0;
          eta_new = (gp[n])[0];
          zeta_new = (gp[n])[1];
          break;
        case 4:
          // face normal in positive zeta direction
          xi_new = (gp[n])[0];
          eta_new = (gp[n])[1];
          zeta_new = +1.0;
          break;
        case 5:
          // face normal in neative xi direction
          xi_new = -1.0;
          eta_new = (gp[n])[1];
          zeta_new = (gp[n])[0];
          break;
        default:
          throw cException(" **** invalid request for face **** ", __FILE__,
                           __LINE__);
          break;
      }

      (gp[n])[0] = xi_new;
      (gp[n])[1] = eta_new;
      (gp[n])[2] = zeta_new;
    }
  } else {
    // ToDo
    // should not get here
    throw cException("Error: Not yet implemented for tetra-elements!", __FILE__,
                     __LINE__);
  }
}

void cElementFluid3d::calculateReactiveSoundPowerSurfaceIntegral(
    Vec fullSolution, int Face, PetscReal &Ps_face_int) {
  // --- Get the shape functions for the face
  cArray3d Nface(getShapeFunctionsFace());
  // --- Get the shape functions for element
  cArray3d N(getShapeFunctions());

  std::vector<short> SurfaceNodes(
      getIndicesOfFaceNodes(Face));  // nodes on one element's face
  std::vector<PetscReal> weights;
  std::vector<cPoint> eps;  // gauss points
  const int nnod_face =
      (int)SurfaceNodes.size();         // number of nodes of one face
  const int nnod = getNumberOfNodes();  // number of nodes of the element
  cElementVector valP_face(nnod_face);
  cElementVector valP(nnod);
  int ngps;
  cVector a1(3);  // vectors describing tangential plane
  cVector a2(3);  // ... at current Gauss point
  cArray3d Nk;    //  shape functions evaluated at nodal points
  PetscInt index;
  PetscReal detJac_face = 0;

  // get the GaussPoints of the face
  getGaussPointsFace(Face, eps, ngps, weights, nnod_face);

  // --- Initialize the shape functions for the velocity calculation
  Nk.initialize(4, nnod, 1);

  // normal vector on face
  cVector normal = getNormalOnFace(Face);

  const PetscReal omega = m_Material->getOmega();  // angle frequency
  const PetscReal rhof = m_Material->getRho();     // fluid density

  // --- extract fluid pressure from solution vector at all the element nodes
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP[k]));
  }

  // --- extract fluid pressure from solution vector at the nodes of the face
  for (int z = 0; z < nnod_face; z++) {
    index = m_Nodes[SurfaceNodes[z]]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP_face[z]));
  }

  for (int n = 0; n < ngps; n++) {
    cElementVector gradPPoint(3);   // grad p at the Gauss Point
    PetscScalar pressurePoint = 0;  // pressure p at the Gauss Point

    a1.setValue(0.0);
    a2.setValue(0.0);

    // --- compute detJac, see e.g.
    //     Boundary Integral Equation Methods for Solids and Fluids
    //     Marc Bonnet
    //     John Wiley & Sons Ltd. (1995)
    //     mech J 111
    /*    for (int z=0; z<nnod_face; z++) {
          for (int d=0; d<3; d++) {
            a1[d] += Nface(1,z,n) * (*m_Nodes[ SurfaceNodes[z] ])[d];
            a2[d] += Nface(2,z,n) * (*m_Nodes[ SurfaceNodes[z] ])[d];
          }
        }
    */
    // by Meike (for curved elements the calculation has to be more detailed)
    for (int k = 0; k < nnod; k++) {
      for (int d = 0; d < 3; d++) {
        if (Face == 3 || Face == 5) {  // d xi
          a1[d] += N(2, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(3, k, n) * (*m_Nodes[k])[d];
        } else if (Face == 1 || Face == 0) {  // d eta
          a1[d] += N(1, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(3, k, n) * (*m_Nodes[k])[d];
        } else if (Face == 4 || Face == 2) {  // d teta
          a1[d] += N(1, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(2, k, n) * (*m_Nodes[k])[d];
        } else
          throw cException("Face nr not known", __FILE__, __LINE__);
      }
    }

    detJac_face = sqrt(a1.abs2() * a2.abs2() - (a1.dot(a2)) * (a1.dot(a2)));

    if (detJac_face <= 0.0) {
      throw cException(
          "detJac_face <= 0.0 in cElementFluid::evaluateSurfaceIntegral()",
          __FILE__, __LINE__);
    }

    // --- Caluculate the gradient of p at the Gauss Point
    computeGradPPoint((eps[n]), Nk, gradPPoint, valP);

    // --- Calculate the pressure value at the Gauss Point
    for (int z = 0; z < nnod_face; z++) {
      pressurePoint += Nface(0, z, n) * valP_face[z];
    }

#ifdef PETSC_USE_COMPLEX

    // --- Caluculate the sound velocity in normal direction at the Gauss Point
    PetscScalar vn = std::complex<PetscReal>(.0, 1.0) /
                     (m_Material->getRhoOmega() * omega) *
                     (gradPPoint[0] * normal[0] + gradPPoint[1] * normal[1] +
                      gradPPoint[2] * normal[2]);

    // --- sound power  (0.5 * because the sound power is related to the
    // effective values of pressure and velocity)
    PetscScalar I_p = 0.5 * pressurePoint * vn;

    // --- abs sound power
    PetscReal I_abs =
        std::sqrt(I_p.real() * I_p.real() + I_p.imag() * I_p.imag());
#else
    // --- Caluculate the sound velocity in normal direction at the Gauss Point
    PetscScalar vn = 1.0 / (m_Material->getRhoOmega() * omega) *
                     (gradPPoint[0] * normal[0] + gradPPoint[1] * normal[1] +
                      gradPPoint[2] * normal[2]);

    // --- sound power  (0.5 * because the sound power is related to the
    // effective values of pressure and velocity)
    PetscScalar I_p = 0.5 * pressurePoint * vn;

    // --- abs sound power
    PetscReal I_abs = std::abs(I_p);
#endif

    Ps_face_int = Ps_face_int + I_abs * weights[n] * detJac_face;
  }
}

void cElementFluid3d::calculateActiveSoundPowerSurfaceIntegral(
    Vec fullSolution, int Face, PetscReal &Ps_face_int) {
  // --- Get the shape functions for the face
  cArray3d Nface(getShapeFunctionsFace());

  // --- Get the shape functions for element
  cArray3d N(getShapeFunctions());

  std::vector<short> SurfaceNodes(
      getIndicesOfFaceNodes(Face));  // nodes on one element's face
  std::vector<PetscReal> weights;
  std::vector<cPoint> eps;  // gauss points
  const int nnod_face =
      (int)SurfaceNodes.size();         // number of nodes of one face
  const int nnod = getNumberOfNodes();  // number of nodes of the element
  cElementVector valP_face(nnod_face);
  cElementVector valP(nnod);
  int ngps;
  cVector a1(3);  // vectors describing tangential plane
  cVector a2(3);  // ... at current Gauss point
  cArray3d Nk;    //  shape functions evaluated at nodal points
  PetscInt index;
  PetscReal detJac_face = 0;
  std::vector<PetscReal> Ps_face_int_sum;

  // get the GaussPoints vor the face
  getGaussPointsFace(Face, eps, ngps, weights, nnod_face);

  // --- Initialize the shape functions for the velocity calculation
  Nk.initialize(4, nnod, 1);

  Ps_face_int_sum.resize(ngps);

  // normal vector on face
  cVector normal = getNormalOnFace(Face);

  const PetscReal omega = m_Material->getOmega();  // angle frequency
  const PetscReal rhof = m_Material->getRho();     // fluid density

  // --- extract fluid pressure from solution vector at all the element nodes
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP[k]));
  }

  // --- extract fluid pressure from solution vector at the nodes of the face
  for (int z = 0; z < nnod_face; z++) {
    index = m_Nodes[SurfaceNodes[z]]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP_face[z]));
  }

  for (int n = 0; n < ngps; n++) {
    cElementVector gradPPoint(3);   // grad p at the Gauss Point
    PetscScalar pressurePoint = 0;  // pressure p at the Gauss Point

    a1.setValue(0.0);
    a2.setValue(0.0);

    // --- compute detJac, see e.g.
    //     Boundary Integral Equation Methods for Solids and Fluids
    //     Marc Bonnet
    //     John Wiley & Sons Ltd. (1995)
    //     mech J 111
    /*    for (int z=0; z<nnod_face; z++) {
          for (int d=0; d<3; d++) {
            a1[d] += Nface(1,z,n) * (*m_Nodes[ SurfaceNodes[z] ])[d];
            a2[d] += Nface(2,z,n) * (*m_Nodes[ SurfaceNodes[z] ])[d];
          }
        }
        PetscReal detJac_face_test = sqrt( a1.abs2() * a2.abs2() -
       (a1.dot(a2))*(a1.dot(a2)));


        a1.setValue(0.0);
        a2.setValue(0.0);
    */
    // by Meike (for curved elements the calculation has to be more detailed)
    for (int k = 0; k < nnod; k++) {
      for (int d = 0; d < 3; d++) {
        if (Face == 3 || Face == 5) {  // d xi
          a1[d] += N(2, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(3, k, n) * (*m_Nodes[k])[d];
        } else if (Face == 1 || Face == 0) {  // d eta
          a1[d] += N(1, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(3, k, n) * (*m_Nodes[k])[d];
        } else if (Face == 4 || Face == 2) {  // d teta
          a1[d] += N(1, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(2, k, n) * (*m_Nodes[k])[d];
        } else
          throw cException("Face nr not known", __FILE__, __LINE__);
      }
    }

    detJac_face = sqrt(a1.abs2() * a2.abs2() - (a1.dot(a2)) * (a1.dot(a2)));

    if (detJac_face <= 0.0) {
      throw cException(
          "detJac_face <= 0.0 in cElementFluid::evaluateSurfaceIntegral()",
          __FILE__, __LINE__);
    }

    // --- Caluculate the gradient of p at the Gauss Point
    computeGradPPoint((eps[n]), Nk, gradPPoint, valP);

    // --- Calculate the pressure value at the Gauss Point
    for (int z = 0; z < nnod_face; z++) {
      pressurePoint += Nface(0, z, n) * valP_face[z];
    }

#ifdef PETSC_USE_COMPLEX

    PetscScalar gradPPoint_n =
        (gradPPoint[0] * normal[0] + gradPPoint[1] * normal[1] +
         gradPPoint[2] * normal[2]);

    // --- Caluculate the sound velocity in normal direction at the Gauss Point
    PetscScalar vn = std::complex<PetscReal>(.0, 1.0) /
                     (m_Material->getRhoOmega() * omega) *
                     (gradPPoint[0] * normal[0] + gradPPoint[1] * normal[1] +
                      gradPPoint[2] * normal[2]);

    // The conjugate complex value of averagevn
    PetscScalar vnconjcomplex =
        vn.real() - std::complex<PetscReal>(.0, 1.0) * vn.imag();

    // --- sound power  (0.5 * because the sound power is related to the
    // effective values of pressure and velocity)
    PetscScalar I_p = 0.5 * pressurePoint * vnconjcomplex;

    // --- abs sound power
    PetscReal I_abs = I_p.real();

#else
    // --- Caluculate the sound velocity in normal direction at the Gauss Point
    PetscScalar vn = 1.0 / (m_Material->getRhoOmega() * omega) *
                     (gradPPoint[0] * normal[0] + gradPPoint[1] * normal[1] +
                      gradPPoint[2] * normal[2]);

    // --- sound power  (0.5 * because the sound power is related to the
    // effective values of pressure and velocity)
    PetscScalar I_p = 0.5 * pressurePoint * vn;

    // --- abs sound power
    PetscReal I_abs = std::abs(I_p);
#endif

    Ps_face_int_sum[n] = I_abs * weights[n] * detJac_face;
  }

  std::sort(Ps_face_int_sum.begin(), Ps_face_int_sum.end(), abs_order);

  for (int i = Ps_face_int_sum.size() - 1; i >= 0; i--)
    Ps_face_int += Ps_face_int_sum[i];
}

void cElementFluid3d::calculateReactiveSoundPowerSurfaceAverage(
    Vec fullSolution, int Face, PetscReal &Ps_face_average) {
  std::vector<short> SurfaceNodes(
      getIndicesOfFaceNodes(Face));  // nodes on one element's face
  const int nnod_face =
      (int)SurfaceNodes.size();         // number of nodes of one face
  const int nnod = getNumberOfNodes();  // number of nodes of the element
  cElementVector valP_face(nnod_face);
  cElementVector valP(nnod);
  cElementVector gradP_face(3 * nnod_face);  // grad p at the face nodes
  cElementVector gradP(3 * nnod);  // grad p at the nodes of the element
  PetscInt index;
  PetscScalar averagePressure = 0;  // average pressure on the face
  PetscScalar averageGradPn =
      0;  // average grad p in nomal direction on the face
  PetscReal area = getFaceArea(Face);

  // normal vector on face
  cVector normal = getNormalOnFace(Face);

  const PetscReal omega = m_Material->getOmega();  // angle frequency
  const PetscReal rhof = m_Material->getRho();     // fluid density

  // --- extract fluid pressure from solution vector at all the element nodes
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP[k]));
  }

  // --- extract fluid pressure from solution vector at the nodes of the face
  for (int z = 0; z < nnod_face; z++) {
    index = m_Nodes[SurfaceNodes[z]]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP_face[z]));

    averagePressure = averagePressure + 1. / nnod_face * valP_face[z];
  }

  // --- calculate the gradient of p at all nodes
  computeGradP(gradP, valP);

  // --- finding gradient of p at the nodes of the face
  for (int z = 0; z < nnod_face; z++) {
    gradP_face[3 * z] = gradP[3 * SurfaceNodes[z]];
    gradP_face[3 * z + 1] = gradP[3 * SurfaceNodes[z] + 1];
    gradP_face[3 * z + 2] = gradP[3 * SurfaceNodes[z] + 2];

    // This is only valid for uncurved elements
    PetscScalar pn = gradP_face[3 * z] * normal[0] +
                     gradP_face[3 * z + 1] * normal[1] +
                     gradP_face[3 * z + 2] * normal[2];
    averageGradPn = averageGradPn + 1. / nnod_face * pn;
  }

#ifdef PETSC_USE_COMPLEX

  // Caluculate the average sound velocity in normal direction
  PetscScalar averagevn = std::complex<PetscReal>(.0, 1.0) /
                          (m_Material->getRhoOmega() * omega) * averageGradPn;

  // average sound power (0.5 * because the sound power is related to the
  // effective values of pressure and velocity)
  PetscScalar I_p = 0.5 * averagePressure * averagevn;

  // abs average sound power
  Ps_face_average =
      std::sqrt(I_p.real() * I_p.real() + I_p.imag() * I_p.imag()) * area;

#else
  // Caluculate the average sound velocity in normal direction
  PetscScalar averagevn =
      1.0 / (m_Material->getRhoOmega() * omega) * averageGradPn;

  // average sound power (0.5 * because the sound power is related to the
  // effective values of pressure and velocity)
  PetscScalar I_p = 0.5 * averagePressure * averagevn;

  // abs average sound power
  Ps_face_average = std::abs(I_p) * area;

#endif
}

void cElementFluid3d::calculateActiveSoundPowerSurfaceAverage(
    Vec fullSolution, int Face, PetscReal &Ps_face_average) {
  std::vector<short> SurfaceNodes(
      getIndicesOfFaceNodes(Face));  // nodes on one element's face
  const int nnod_face =
      (int)SurfaceNodes.size();         // number of nodes of one face
  const int nnod = getNumberOfNodes();  // number of nodes of the element
  cElementVector valP_face(nnod_face);
  cElementVector valP(nnod);
  cElementVector gradP_face(3 * nnod_face);  // grad p at the face nodes
  cElementVector gradP(3 * nnod);  // grad p at the nodes of the element
  PetscInt index;
  PetscScalar averagePressure = 0;  // average pressure on the face
  PetscScalar averageGradPn =
      0;  // average grad p in nomal direction on the face
  PetscReal area = getFaceArea(Face);

  // normal vector on face
  cVector normal = getNormalOnFace(Face);

  const PetscReal omega = m_Material->getOmega();  // angle frequency
  const PetscReal rhof = m_Material->getRho();     // fluid density

  // --- extract fluid pressure from solution vector at all the element nodes
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP[k]));
  }

  // --- extract fluid pressure from solution vector at the nodes of the face
  for (int z = 0; z < nnod_face; z++) {
    index = m_Nodes[SurfaceNodes[z]]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP_face[z]));

    averagePressure = averagePressure + 1. / nnod_face * valP_face[z];
  }

  // --- calculate the gradient of p at all nodes
  computeGradP(gradP, valP);

  // --- finding gradient of p at the nodes of the face
  for (int z = 0; z < nnod_face; z++) {
    gradP_face[3 * z] = gradP[3 * SurfaceNodes[z]];
    gradP_face[3 * z + 1] = gradP[3 * SurfaceNodes[z] + 1];
    gradP_face[3 * z + 2] = gradP[3 * SurfaceNodes[z] + 2];

    // This is only valid for uncurved elements
    PetscScalar pn = gradP_face[3 * z] * normal[0] +
                     gradP_face[3 * z + 1] * normal[1] +
                     gradP_face[3 * z + 2] * normal[2];
    averageGradPn = averageGradPn + 1. / nnod_face * pn;
  }

#ifdef PETSC_USE_COMPLEX

  // Caluculate the average sound velocity in normal direction
  PetscScalar averagevn = std::complex<PetscReal>(.0, 1.0) /
                          (m_Material->getRhoOmega() * omega) * averageGradPn;

  // The conjugate complex value of averagevn
  PetscScalar vnconjcomplex =
      averagevn.real() - std::complex<PetscReal>(.0, 1.0) * averagevn.imag();

  // average sound power (0.5 * because the sound power is related to the
  // effective values of pressure and velocity)
  PetscScalar I_p = 0.5 * averagePressure * vnconjcomplex;

  // abs average sound power
  Ps_face_average = I_p.real() * area;

#else
  // Caluculate the average sound velocity in normal direction
  PetscScalar averagevn =
      1.0 / (m_Material->getRhoOmega() * omega) * averageGradPn;

  // average sound power (0.5 * because the sound power is related to the
  // effective values of pressure and velocity)
  PetscScalar I_p = 0.5 * averagePressure * averagevn;

  // abs average sound power
  Ps_face_average = std::abs(I_p) * area;

#endif
}

void cElementFluid3d::calculateActiveSoundPowerSurfaceIntegralGivenVn(
    Vec fullSolution, int Face, PetscReal &Ps_face_int, PetscReal GivenVn) {
  // --- Get the shape functions for the face
  cArray3d Nface(getShapeFunctionsFace());

  // --- Get the shape functions for element
  cArray3d N(getShapeFunctions());

  std::vector<short> SurfaceNodes(
      getIndicesOfFaceNodes(Face));  // nodes on one element's face
  std::vector<PetscReal> weights;
  std::vector<cPoint> eps;  // gauss points
  const int nnod_face =
      (int)SurfaceNodes.size();         // number of nodes of one face
  const int nnod = getNumberOfNodes();  // number of nodes of the element
  cElementVector valP_face(nnod_face);
  cElementVector valP(nnod);
  int ngps;
  cVector a1(3);  // vectors describing tangential plane
  cVector a2(3);  // ... at current Gauss point
  cArray3d Nk;    //  shape functions evaluated at nodal points
  PetscInt index;
  PetscReal detJac_face = 0;
  std::vector<PetscReal> Ps_face_int_sum;

  // get the GaussPoints vor the face
  getGaussPointsFace(Face, eps, ngps, weights, nnod_face);

  // --- Initialize the shape functions for the velocity calculation
  Nk.initialize(4, nnod, 1);

  Ps_face_int_sum.resize(ngps);

  // normal vector on face
  cVector normal = getNormalOnFace(Face);

  const PetscReal omega = m_Material->getOmega();  // angle frequency
  const PetscReal rhof = m_Material->getRho();     // fluid density

  // --- extract fluid pressure from solution vector at all the element nodes
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP[k]));
  }

  // --- extract fluid pressure from solution vector at the nodes of the face
  for (int z = 0; z < nnod_face; z++) {
    index = m_Nodes[SurfaceNodes[z]]->getGlobalRow(fluid);
    VecGetValues(fullSolution, 1, &index, &(valP_face[z]));
  }

  for (int n = 0; n < ngps; n++) {
    PetscScalar pressurePoint = 0;  // pressure p at the Gauss Point

    a1.setValue(0.0);
    a2.setValue(0.0);

    // --- compute detJac, see e.g.
    //     Boundary Integral Equation Methods for Solids and Fluids
    //     Marc Bonnet
    //     John Wiley & Sons Ltd. (1995)
    //     mech J 111
    /*    for (int z=0; z<nnod_face; z++) {
          for (int d=0; d<3; d++) {
            a1[d] += Nface(1,z,n) * (*m_Nodes[ SurfaceNodes[z] ])[d];
            a2[d] += Nface(2,z,n) * (*m_Nodes[ SurfaceNodes[z] ])[d];
          }
        }
        PetscReal detJac_face_test = sqrt( a1.abs2() * a2.abs2() -
       (a1.dot(a2))*(a1.dot(a2)));


        a1.setValue(0.0);
        a2.setValue(0.0);
    */
    // by Meike (for curved elements the calculation has to be more detailed)
    for (int k = 0; k < nnod; k++) {
      for (int d = 0; d < 3; d++) {
        if (Face == 3 || Face == 5) {  // d xi
          a1[d] += N(2, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(3, k, n) * (*m_Nodes[k])[d];
        } else if (Face == 1 || Face == 0) {  // d eta
          a1[d] += N(1, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(3, k, n) * (*m_Nodes[k])[d];
        } else if (Face == 4 || Face == 2) {  // d teta
          a1[d] += N(1, k, n) * (*m_Nodes[k])[d];
          a2[d] += N(2, k, n) * (*m_Nodes[k])[d];
        } else
          throw cException("Face nr not known", __FILE__, __LINE__);
      }
    }
    detJac_face = sqrt(a1.abs2() * a2.abs2() - (a1.dot(a2)) * (a1.dot(a2)));

    if (detJac_face <= 0.0) {
      throw cException(
          "detJac_face <= 0.0 in cElementFluid::evaluateSurfaceIntegral()",
          __FILE__, __LINE__);
    }

    // --- Calculate the pressure value at the Gauss Point
    for (int z = 0; z < nnod_face; z++) {
      pressurePoint += Nface(0, z, n) * valP_face[z];
    }

#ifdef PETSC_USE_COMPLEX

    // PetscScalar gradPPoint_n = (gradPPoint[0]*normal[0] +
    // gradPPoint[1]*normal[1] + gradPPoint[2]*normal[2]);

    // --- Caluculate the sound velocity in normal direction at the Gauss Point
    PetscScalar vn = GivenVn;

    // The conjugate complex value of averagevn
    PetscScalar vnconjcomplex =
        vn.real() - std::complex<PetscReal>(.0, 1.0) * vn.imag();

    // --- sound power  (0.5 * because the sound power is related to the
    // effective values of pressure and velocity)
    PetscScalar I_p = 0.5 * pressurePoint * vnconjcomplex;

    // --- abs sound power
    PetscReal I_abs = I_p.real();

#else
    // --- Caluculate the sound velocity in normal direction at the Gauss Point
    PetscScalar vn =
        GivenVn;  // 1.0 / (m_Material->getRhoOmega() * omega) *
                  // (gradPPoint[0]*normal[0] + gradPPoint[1]*normal[1] +
                  // gradPPoint[2]*normal[2]);

    // --- sound power  (0.5 * because the sound power is related to the
    // effective values of pressure and velocity)
    PetscScalar I_p = 0.5 * pressurePoint * vn;

    // --- abs sound power
    PetscReal I_abs = std::abs(I_p);
#endif

    Ps_face_int_sum[n] = I_abs * weights[n] * detJac_face;
  }

  std::sort(Ps_face_int_sum.begin(), Ps_face_int_sum.end(), abs_order);

  for (int i = Ps_face_int_sum.size() - 1; i >= 0; i--)
    Ps_face_int += Ps_face_int_sum[i];
}

PetscReal cElementFluid3d::getFaceArea(const short &face) {
  // --- Get the shape functions for the face
  cArray3d Nface(getShapeFunctionsFace());
  std::vector<short> SurfaceNodes(
      getIndicesOfFaceNodes(face));  // nodes on one element's face
  const int nnod_face =
      (int)SurfaceNodes.size();  // number of nodes of one face
  int level;                     // number of Gauss points to use
  int ngps;
  cVector a1(3);  // vectors describing tangential plane
  cVector a2(3);  // ... at current Gauss point
  std::vector<PetscReal> weights;

  PetscReal area = 0;
  PetscReal detJac_face = 0;

  if (nnod_face == 3) {
    level = 1;
    ngps = level;
    weights.resize(ngps);
    for (int n = 0; n < ngps; n++) {
      weights[n] = 0.5 * m_GaussPoints.getGaussWeightTria(level, n);
    }
  } else if (nnod_face == 4) {
    level = 2;
    ngps = level * level;
    weights.resize(ngps);
    for (int n = 0; n < ngps; n++) {
      weights[n] = m_GaussPoints.getGaussWeight2D(level, n);
    }
  } else if (nnod_face == 9) {
    level = 3;
    ngps = level * level;
    weights.resize(ngps);
    for (int n = 0; n < ngps; n++) {
      weights[n] = m_GaussPoints.getGaussWeight2D(level, n);
    }
  } else {
    // should not get here
    throw cException(
        "Error: do not know what to do for this value of nnod_face !", __FILE__,
        __LINE__);
  }

  for (int n = 0; n < ngps; n++) {
    a1.setValue(0.0);
    a2.setValue(0.0);

    // --- compute detJac, see e.g.
    //     Boundary Integral Equation Methods for Solids and Fluids
    //     Marc Bonnet
    //     John Wiley & Sons Ltd. (1995)
    //     mech J 111
    for (int z = 0; z < nnod_face; z++) {
      for (int d = 0; d < 3; d++) {
        a1[d] += Nface(1, z, n) * (*m_Nodes[SurfaceNodes[z]])[d];
        a2[d] += Nface(2, z, n) * (*m_Nodes[SurfaceNodes[z]])[d];
      }
    }

    detJac_face = sqrt(a1.abs2() * a2.abs2() - (a1.dot(a2)) * (a1.dot(a2)));

    area = area + weights[n] * detJac_face;
  }
  return area;
}

// calculates phase shift of pressure and velocity
void cElementFluid3d::computePhaseShiftNodes(Vec FullSolution,
                                             Vec &phaseShiftAverageVec,
                                             Vec &phaseShiftVec) {
  const int nnod = getNumberOfNodes();

  const PetscReal omega = m_Material->getOmega();  // angle frequency
  const PetscReal rhof =
      m_Material->getRho();  // frequency depdendent fluid density
  //  const PetscReal cf    = m_Material->getCf();		 // frequency depdendent
  //  speed of sound

  cElementVector gradP(3 * nnod);
  cElementVector phaseShiftAverage(nnod);
  cElementVector phaseShift(3 * nnod);
  cElementVector valP(nnod);
  PetscInt index, seqId;
  std::vector<PetscScalar> v(3);  // velocitys in x,y,z direction
  std::vector<PetscReal> angle_v(3);

  PetscInt *indicesPSA = new PetscInt[nnod];
  PetscInt *indicesPS = new PetscInt[3 * nnod];

  // --- extract fluid pressure from solution vector
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    seqId = m_Nodes[k]->getGlobalSeqId();

    indicesPSA[k] = index;
    indicesPS[3 * k] = cstNumberOfStressDofs * seqId + gradpx;
    indicesPS[3 * k + 1] = cstNumberOfStressDofs * seqId + gradpy;
    indicesPS[3 * k + 2] = cstNumberOfStressDofs * seqId + gradpz;

    VecGetValues(FullSolution, 1, &index, &(valP[k]));
  }

  computeGradP(gradP, valP);

  // --- calculation of PhaseShift
  for (int k = 0; k < nnod; k++) {
    // divide through nrValuesAtNode to get the average value and not the summ
    // at the nodes
    PetscInt nrValuesAtNode = m_Nodes[k]->getnrAdjacentElements(gradpx);

    if (nrValuesAtNode < 1)
      trace(
          "ERROR: nrValuesAtNode can not be zero, check the use of "
          "cOutput::countElementsAtNodes ... !!!!!!!!!!");

#ifdef PETSC_USE_COMPLEX

    // velocitys in x,y,z direction
    for (int i = 0; i < 3; i++)
      v[i] =
          std::complex<PetscReal>(.0, 1.0) / (rhof * omega) * gradP[3 * k + i];

    for (int i = 0; i < 3; i++) {
      PetscReal rad_v = std::sqrt((v[i]).real() * (v[i]).real() +
                                  (v[i]).imag() * (v[i]).imag());

      if (rad_v > 0.)
        angle_v[i] = std::acos((v[i]).real() / rad_v);
      else
        angle_v[i] = 0.;

      if ((v[i]).imag() < 0.) angle_v[i] = 2 * M_PI - angle_v[i];
    }

    PetscReal angle_p = 0.;
    PetscReal angle_tmp = 0.;
    PetscReal rad_p = std::sqrt((valP[k]).real() * (valP[k]).real() +
                                (valP[k]).imag() * (valP[k]).imag());
    if (rad_p > 0.)
      angle_p = std::acos((valP[k]).real() / rad_p);
    else
      angle_p = 0.;

    if ((valP[k]).imag() < 0.) angle_p = 2 * M_PI - angle_p;

    // divide through nrValuesAtNode to get the average value and not the summ
    // at the nodes
    for (int i = 0; i < 3; i++) {
      angle_tmp = std::abs(angle_p - angle_v[i]);
      if (angle_tmp > M_PI) angle_tmp = 2 * M_PI - angle_tmp;

      phaseShift[3 * k + i] = 1. / nrValuesAtNode * angle_tmp;

      phaseShiftAverage[k] =
          phaseShiftAverage[k] + 1. / (3. * nrValuesAtNode) * angle_tmp;
    }

#else
    phaseShift[3 * k] = 0;
    phaseShift[3 * k + 1] = 0;
    phaseShift[3 * k + 2] = 0;

    trace(
        "ERROR: phase shift output in frequency range without complex values "
        "makes no sense ... !!!!!!!!!!");
#endif
  }

  // --- add to global vectors
  VecSetValues(phaseShiftAverageVec, nnod, indicesPSA, &(phaseShiftAverage[0]),
               ADD_VALUES);
  VecSetValues(phaseShiftVec, 3 * nnod, indicesPS, &(phaseShift[0]),
               ADD_VALUES);

  // --- delete the arrays
  delete[] indicesPS;
  delete[] indicesPSA;
}

void cElementFluid3d::computeEDandIntNodes(Vec FullSolution, Vec &EnergyDensity,
                                           Vec &intensity) {
  const int nnod = getNumberOfNodes();

  const PetscReal omega = m_Material->getOmega();  // angle frequency
  const PetscReal rhof =
      m_Material->getRho();  // frequency depdendent fluid density
  const PetscReal cf =
      m_Material->getCf();  // frequency depdendent speed of sound

  cElementVector gradP(3 * nnod);
  cElementVector en_dens(nnod);
  cElementVector intens(3 * nnod);
  cElementVector valP(nnod);
  PetscInt index, seqId;
  std::vector<PetscScalar> v(3);  // velocitys in x,y,z direction
  PetscScalar v_mag = 0;          // magnitude of the velocity
  PetscReal p2_abs = 0;
  PetscReal v_mag2_abs = 0;

  PetscInt *indicesED = new PetscInt[nnod];
  PetscInt *indicesI = new PetscInt[3 * nnod];

  // --- extract fluid pressure from solution vector
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    seqId = m_Nodes[k]->getGlobalSeqId();

    indicesED[k] = index;
    indicesI[3 * k] = cstNumberOfStressDofs * seqId + gradpx;
    indicesI[3 * k + 1] = cstNumberOfStressDofs * seqId + gradpy;
    indicesI[3 * k + 2] = cstNumberOfStressDofs * seqId + gradpz;

    VecGetValues(FullSolution, 1, &index, &(valP[k]));
  }

  computeGradP(gradP, valP);

  // --- calculation of EnergyDensity
  for (int k = 0; k < nnod; k++) {
    // divide through nrValuesAtNode to get the average value and not the summ
    // at the nodes
    PetscInt nrValuesAtNode = m_Nodes[k]->getnrAdjacentElements(gradpx);

    if (nrValuesAtNode < 1)
      trace(
          "ERROR: nrValuesAtNode can not be zero, check the use of "
          "cOutput::countElementsAtNodes ... !!!!!!!!!!");

#ifdef PETSC_USE_COMPLEX

    // velocitys in x,y,z direction
    for (int i = 0; i < 3; i++)
      v[i] =
          std::complex<PetscReal>(.0, 1.0) / (rhof * omega) * gradP[3 * k + i];

    // magnitude of the velocity
    v_mag = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

    p2_abs = (valP[k]).real() * (valP[k]).real() +
             (valP[k]).imag() * (valP[k]).imag();
    v_mag2_abs = v_mag.real() * v_mag.real() + v_mag.imag() * v_mag.imag();

    PetscScalar I_0 =
        0.5 * valP[k] *
        ((v[0]).real() - std::complex<PetscReal>(.0, 1.0) * (v[0]).imag());
    PetscScalar I_1 =
        0.5 * valP[k] *
        ((v[1]).real() - std::complex<PetscReal>(.0, 1.0) * (v[1]).imag());
    PetscScalar I_2 =
        0.5 * valP[k] *
        ((v[2]).real() - std::complex<PetscReal>(.0, 1.0) * (v[2]).imag());

    // divide through nrValuesAtNode to get the average value and not the summ
    // at the nodes
    intens[3 * k] = 1. / nrValuesAtNode * I_0.real();
    intens[3 * k + 1] = 1. / nrValuesAtNode * I_1.real();
    intens[3 * k + 2] = 1. / nrValuesAtNode * I_2.real();

#else
    // velocitys in x,y,z direction
    v[0] = 1.0 / (rhof * omega) * gradP[3 * k];
    v[1] = 1.0 / (rhof * omega) * gradP[3 * k + 1];
    v[2] = 1.0 / (rhof * omega) * gradP[3 * k + 2];

    // magnitude of the velocity
    v_mag = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

    p2_abs = valP[k] * valP[k];
    v_mag2_abs = v_mag * v_mag;

    intens[3 * k] = 0;
    intens[3 * k + 1] = 0;
    intens[3 * k + 2] = 0;

    trace(
        "ERROR: insensity output in frequency range without complex values "
        "makes no sense ... !!!!!!!!!!");

#endif

    // energy density at the node k
    // divide through nrValuesAtNode to get the average value and not the summ
    // at the nodes
    en_dens[k] = 1. / nrValuesAtNode *
                 (p2_abs / (4 * rhof * cf * cf) + rhof * v_mag2_abs / 4);
  }

  // --- add to global vectors
  VecSetValues(EnergyDensity, nnod, indicesED, &(en_dens[0]), ADD_VALUES);

  VecSetValues(intensity, 3 * nnod, indicesI, &(intens[0]), ADD_VALUES);

  // --- delete the arrays
  delete[] indicesED;
  delete[] indicesI;
}
