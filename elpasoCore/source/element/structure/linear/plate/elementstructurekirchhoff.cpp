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

#include "elementstructurekirchhoff.h"


cElementStructureKirchhoff::cElementStructureKirchhoff(short NumberOfNodes, short NumberOfGaussPoints) :
  cElementStructureLinear(NumberOfNodes, 4, NumberOfGaussPoints)
{
  // leer
}

cElementStructureKirchhoff::cElementStructureKirchhoff(const cElementStructureKirchhoff &other) :
  cElementStructureLinear(other)
{
  // leer
}


cElementStructureKirchhoff::~cElementStructureKirchhoff()
{
  // leer
}


std::vector<eKnownDofs> cElementStructureKirchhoff::getDofs(void) const
{
  std::vector<eKnownDofs> res(4);

  res[0] = disp_x3;
  res[1] = disp_dwdx;
  res[2] = disp_dwdy;
  res[3] = disp_dwdxy;

  return res;
}

void cElementStructureKirchhoff::computeStressesCells(Vec &fullSolution, Vec &stresses, Vec &stressesSec2)
{
  // Stress calculation for Kirchhoff plate [Uni Siegen VL Baustatik Kap. 8.3]

  const int nnod  = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();
  int nn=4; // stress calculation only at 4 corner nodes
  cElementVector localSolution(nnod*ndofs);
  getLocalSolutionFromFullSolution(&fullSolution,localSolution);

  cArray3d  Nk;           //  shape functions evaluated at nodal points
  Nk.initialize(3, nnod, 1);

  cMatrix Bm(3, nn*3);       // Matrix of Ansatzfunctions

  cElementVector localSolutionKirch4(nn*ndofs); // reduce nodes to 4
  cElementVector relevantLocalSolutionKirch4(nn*3); // reduce nodes to 4

  for(int j=0;j<nn*ndofs;j++) localSolutionKirch4[j] = localSolution[j]; // copying the local solution from Kirch4 to Kirch4 (maybe in the future, higher order follows..)

  relevantLocalSolutionKirch4[0] = localSolutionKirch4[1]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionKirch4[1] = localSolutionKirch4[2];
  relevantLocalSolutionKirch4[2] = 0.;//localSolutionKirch4[3];
  relevantLocalSolutionKirch4[3] = localSolutionKirch4[5]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionKirch4[4] = localSolutionKirch4[6];
  relevantLocalSolutionKirch4[5] = 0.;//localSolutionKirch4[7];
  relevantLocalSolutionKirch4[6] = localSolutionKirch4[9]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionKirch4[7] = localSolutionKirch4[10];
  relevantLocalSolutionKirch4[8] = 0.;//localSolutionKirch4[11];
  relevantLocalSolutionKirch4[9] = localSolutionKirch4[13]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionKirch4[10] = localSolutionKirch4[14];
  relevantLocalSolutionKirch4[11] = 0.;//localSolutionKirch4[15];

  //cElementVector strain(6);   //local strain
  cElementVector stress(6);   //local stresses
  cElementVector strain2d(3); //relevant local strain
  cElementVector stress2d(3); //relevant local strain

  PetscReal Nx; // N,x
  PetscReal Ny; // N,y
  cPoint    ep;

  // -------------------------------------------------------------------------
  //  Elastizitaetsmatrix bestimmen
  // -------------------------------------------------------------------------
  if(m_history!=NULL) m_Material->setHistory(m_history);
  cElementMatrix Cmps(3,3);
  m_Material->setupCps(Cmps);

  // -- get single node of element. The stress values are accurate at centroid of Kirchhoff plate
  ep[0] = 0.;
  ep[1] = 0.;

  for (int k=0; k<nn; k++)
  {
    // -- evaluate shape functions
    Nk(0,k,0) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nn, k, N_fun, ep);
    Nk(1,k,0) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nn, k, N_xi, ep);
    Nk(2,k,0) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nn, k, N_eta, ep);
  }

  // -----------------------------------------------------------------------
  //  Aufstellen der Jacobi-Matrix
  // -----------------------------------------------------------------------
  setupJacobian2D(Nk, 0); //has to be 0, due to the reduced approach
  invertJacobian2D();

  // -----------------------------------------------------------------------
  //  Aufstellen der Matrix der Ansatzfkt.
  // -----------------------------------------------------------------------
  for (int k=0;k<nn;k++)
  {
    Nx = (invJac(0,0)*Nk(1,k,0) + invJac(0,1)*Nk(2,k,0)); // dN/dx
    Ny = (invJac(1,0)*Nk(1,k,0) + invJac(1,1)*Nk(2,k,0)); // dN/dy

    Bm(0, 0+k*3) = Nx;
    Bm(1, 1+k*3) = Ny;

    Bm(2, 0+k*3) = Ny;
    Bm(2, 1+k*3) = Nx;
  }

  //get local strain
  infam::mult(Bm,relevantLocalSolutionKirch4,strain2d);
  const PetscReal thickness = m_Material->getT(this);
  strain2d[0] = -0.5*thickness*strain2d[0];
  strain2d[1] = -0.5*thickness*strain2d[1];
  strain2d[2] = -0.5*thickness*strain2d[2];

  //get local stress
  infam::mult(Cmps,strain2d,stress2d);

  stress[0] = stress2d[0];
  stress[1] = stress2d[1];
  stress[2] = 0.;
  stress[3] = stress2d[2];
  stress[4] = 0.;
  stress[5] = 0.;

  if(m_history!=NULL)
  {
    //if element has yielding get absolute stress from element history
    if (m_history->getYielding() == true && m_history->hasElasticValues() == true)
    {
      //std::cout<<"Element "<<getId()<<" getAbsoluteStress(stress)\n";
      m_history->getAbsoluteStress(stress);
      //std::cout<<"stress \n"<<stress<<"\n";
    }
  }

  //add element stresses to global stress vector
  PetscInt pos;

  for(int i=0; i<6; ++i)
  {
    pos=(this->getId0n())*9+i;
    VecSetValue(stresses,pos,stress[i],INSERT_VALUES);
  }
}

// --- compute plate velocitys at nodes
void cElementStructureKirchhoff::computeVi(PetscReal* n_vis, std::vector<PetscReal>* vis, const Vec &Solution)
{  
  if (m_ElementLoads.size() == 0) {
	trace( "********************************************************" );
	trace( "* computeVi can only be computed with existing ElementLoads *" );
	trace( "********************************************************" );
    return;
  }
//  if (this->getId() != 66)
//	  return;

  const int      nnod = getNumberOfNodes();
  PetscInt       index;
  cElementVector Us(3 * nnod);
  const PetscReal   omega = m_Material->getOmega();


  // --- get displacements
  for (int k=0; k<nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(disp_x1);
    VecGetValues(Solution, 1, &index, &(Us[3*k+0]));

    index = m_Nodes[k]->getGlobalRow(disp_x2);
    VecGetValues(Solution, 1, &index, &(Us[3*k+1]));

    index = m_Nodes[k]->getGlobalRow(disp_x3);
    VecGetValues(Solution, 1, &index, &(Us[3*k+2]));
  }


  std::multimap<short, cElementLoad *>::const_iterator itLoads;
  for (itLoads=m_ElementLoads.begin(); itLoads!=m_ElementLoads.end(); itLoads++)
  {
    const int nnod = this->getNumberOfNodes();

	// --- get Displacements for the considered faces
    for (int k=0; k<nnod; k++) {
	  PetscScalar ui_val_x = Us[3*k+0];
	  PetscScalar ui_val_y = Us[3*k+1];
	  PetscScalar ui_val_z = Us[3*k+2];

// Meike test
//	  if (m_Nodes[k]->getId() == 72) {

	   // --- compute the velocitys
#ifdef PETSC_USE_COMPLEX
	  (*vis)[0] += omega * std::sqrt(ui_val_x.real()*ui_val_x.real() + ui_val_x.imag()*ui_val_x.imag());
	  (*vis)[1] += omega * std::sqrt(ui_val_y.real()*ui_val_y.real() + ui_val_y.imag()*ui_val_y.imag());
	  (*vis)[2] += omega * std::sqrt(ui_val_z.real()*ui_val_z.real() + ui_val_z.imag()*ui_val_z.imag());
// get the displacements
/*	  (*vis)[0] += std::sqrt(ui_val_x.real()*ui_val_x.real() + ui_val_x.imag()*ui_val_x.imag());
	  (*vis)[1] += std::sqrt(ui_val_y.real()*ui_val_y.real() + ui_val_y.imag()*ui_val_y.imag());
	  (*vis)[2] += std::sqrt(ui_val_z.real()*ui_val_z.real() + ui_val_z.imag()*ui_val_z.imag());*/

#else
	  (*vis)[0] += omega * ui_val_x;
	  (*vis)[1] += omega * ui_val_y;
	  (*vis)[2] += omega * ui_val_z;

/*	  (*vis)[0] +=  ui_val_x;
	  (*vis)[1] +=  ui_val_y;
	  (*vis)[2] +=  ui_val_z;*/
#endif
      *n_vis += 1.;
//	  }
    }
  }
}
