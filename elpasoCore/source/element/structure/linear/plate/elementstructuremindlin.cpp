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

#include "elementstructuremindlin.h"

cElementStructureMindlin::cElementStructureMindlin(short NumberOfNodes, short NumberOfGaussPoints) :
cElementStructureLinear(NumberOfNodes, 3, NumberOfGaussPoints)
{
  // empty
}

cElementStructureMindlin::cElementStructureMindlin(const cElementStructureMindlin &other) :
cElementStructureLinear(other)
{
  // empty
}

cElementStructureMindlin::~cElementStructureMindlin()
{
  // empty
}


std::vector<eKnownDofs> cElementStructureMindlin::getDofs(void) const
{
  std::vector<eKnownDofs> res(3);

  res[0] = disp_x3;
  res[1] = disp_w1;
  res[2] = disp_w2;

  return res;
}

void cElementStructureMindlin::computeStressesCells(Vec &fullSolution, Vec &stresses, Vec &stressesSec2)
{
  // Stress calculation for Mindlin plate [Uni Siegen VL Baustatik Kap. 8.3]

  const int nnod  = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();
  int nn=4; // stress calculation only at 4 corner nodes
  cElementVector localSolution(nnod*ndofs);
  getLocalSolutionFromFullSolution(&fullSolution,localSolution);

  cArray3d  Nk;           //  shape functions evaluated at nodal points
  Nk.initialize(3, nnod, 1);

  cMatrix Bm(3, nn*3);       // Matrix of Ansatzfunctions

  cElementVector localSolutionDSG4(nn*ndofs); // reduce nodes to 4
  cElementVector relevantLocalSolutionDSG4(nn*3); // reduce nodes to 4

  for(int j=0;j<nn*ndofs;j++) localSolutionDSG4[j] = localSolution[j]; // copying the local solution from DSG4 to DSG4 (maybe in the future, higher order follows..)

  relevantLocalSolutionDSG4[0] = localSolutionDSG4[2]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionDSG4[1] = localSolutionDSG4[1];
  relevantLocalSolutionDSG4[2] = localSolutionDSG4[2];
  relevantLocalSolutionDSG4[3] = localSolutionDSG4[5]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionDSG4[4] = localSolutionDSG4[4];
  relevantLocalSolutionDSG4[5] = localSolutionDSG4[5];
  relevantLocalSolutionDSG4[6] = localSolutionDSG4[8]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionDSG4[7] = localSolutionDSG4[7];
  relevantLocalSolutionDSG4[8] = localSolutionDSG4[8];
  relevantLocalSolutionDSG4[9] = localSolutionDSG4[11]; // Only direct given derivatives necessary , dof w not necessary
  relevantLocalSolutionDSG4[10] = localSolutionDSG4[10];
  relevantLocalSolutionDSG4[11] = localSolutionDSG4[11];

  cElementVector stress(6);   //local stresses
  cElementVector stressSec2(6);   //local stresses - second section (lower side)
  cElementVector strain2d(3); //relevant local strain
  cElementVector strain2dSec2(3); //relevant local strain - second section (lower side)
  cElementVector stress2d(3); //relevant local stress
  cElementVector stress2dSec2(3); //relevant local stress - second section (lower side)

  PetscReal Nx; // N,x
  PetscReal Ny; // N,y
  cPoint    ep;

  // -------------------------------------------------------------------------
  //  Elastizitaetsmatrix bestimmen
  // -------------------------------------------------------------------------
  if(m_history!=NULL) m_Material->setHistory(m_history);
  cElementMatrix Cm(3,3);
  m_Material->setupCm(Cm);

  // -- get single node of element. The stress values are accurate at centroid of Mindlin plate
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
    Bm(2, 2+k*3) = Ny;
  }

  //get local strain
  infam::mult(Bm,relevantLocalSolutionDSG4,strain2d);

  const PetscReal thickness = m_Material->getT(this);
  strain2dSec2[0] = -0.5*thickness*strain2d[0];
  strain2dSec2[1] = +0.5*thickness*strain2d[1];
  strain2dSec2[2] = -0.5*2.0*thickness*strain2d[2];

  strain2d[0] = 0.5*thickness*strain2d[0];
  strain2d[1] = -0.5*thickness*strain2d[1];
  strain2d[2] = 0.5*2.0*thickness*strain2d[2];

  //get local stress
  infam::mult(Cm,strain2d,stress2d);
  infam::mult(Cm,strain2dSec2,stress2dSec2);

  stress[0] = stress2d[0];
  stress[1] = stress2d[1];
  stress[2] = 0.;
  stress[3] = stress2d[2];
  stress[4] = 0.;
  stress[5] = 0.;

  stressSec2[0] = stress2dSec2[0];
  stressSec2[1] = stress2dSec2[1];
  stressSec2[2] = 0.;
  stressSec2[3] = stress2dSec2[2];
  stressSec2[4] = 0.;
  stressSec2[5] = 0.;

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
    VecSetValue(stressesSec2,pos,stressSec2[i],INSERT_VALUES);
  }
}

void cElementStructureMindlin::getStressesCells(cElementVector &relevantLocalSolutionDSG4, cElementVector &stress2d, cElementVector &stress2dSec2)
{
  // Stress calculation for Mindlin plate [Uni Siegen VL Baustatik Kap. 8.3]

  const int nnod  = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();
  int nn=4; // stress calculation only at 4 corner nodes

  cArray3d  Nk;           //  shape functions evaluated at nodal points
  Nk.initialize(3, nnod, 1);

  cMatrix Bm(3, nn*3);       // Matrix of Ansatzfunctions

  cElementVector stress(6);   //local stresses
  //cElementVector stressSec2(6);   //local stresses - second section (lower side)
  cElementVector strain2d(3); //relevant local strain
  cElementVector strain2dSec2(3); //relevant local strain - second section (lower side)
  //cElementVector stress2d(3); //relevant local strain
  //cElementVector stress2dSec2(3); //relevant local strain - second section (lower side)

  PetscReal Nx; // N,x
  PetscReal Ny; // N,y
  cPoint    ep;

  // -------------------------------------------------------------------------
  //  Elastizitaetsmatrix bestimmen
  // -------------------------------------------------------------------------
  if(m_history!=NULL) m_Material->setHistory(m_history);
  cElementMatrix Cm(3,3);
  m_Material->setupCm(Cm);

  // -- get single node of element. The stress values are accurate at centroid of Mindlin plate
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
    Bm(2, 2+k*3) = Ny;
  }

  //get local strain
  infam::mult(Bm,relevantLocalSolutionDSG4,strain2d);

  const PetscReal thickness = m_Material->getT(this);
  strain2dSec2[0] = -0.5*thickness*strain2d[0];
  strain2dSec2[1] = +0.5*thickness*strain2d[1];
  strain2dSec2[2] = -0.5*2.0*thickness*strain2d[2];

  strain2d[0] = 0.5*thickness*strain2d[0];
  strain2d[1] = -0.5*thickness*strain2d[1];
  strain2d[2] = 0.5*2.0*thickness*strain2d[2];

  //get local stress
  infam::mult(Cm,strain2d,stress2d);
  infam::mult(Cm,strain2dSec2,stress2dSec2);

  /*stress[0] = stress2d[0];
  stress[1] = stress2d[1];
  stress[2] = 0.;
  stress[3] = stress2d[2];
  stress[4] = 0.;
  stress[5] = 0.;

  stressSec2[0] = stress2dSec2[0];
  stressSec2[1] = stress2dSec2[1];
  stressSec2[2] = 0.;
  stressSec2[3] = stress2dSec2[2];
  stressSec2[4] = 0.;
  stressSec2[5] = 0.;

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
    VecSetValue(stressesSec2,pos,stressSec2[i],INSERT_VALUES);
  }*/
}

void cElementStructureMindlin::assembleMassMatrix(cElementMatrix &MM)
{
  PetscReal        wdJ;       // Gewicht * detJ
  const short ngp = getNumberOfGaussPoints();
  cArray3d    N( getShapeFunctions() );


  // -------------------------------------------------------------------------
  //  benoetigte Materialparameter
  // -------------------------------------------------------------------------
  const int  nnod  = getNumberOfNodes();  // Anzahl der Knoten im Element
  const int  ndofs = 3;                   // Anzahl Freiheitsgrade je Knoten

  const PetscReal rho = m_Material->getRho();  // Dichte des Materials
  const PetscReal h   = m_Material->getT(this);    // Dicke des Plattenelements
  const PetscReal ftm = h*h*h/12.0;
  const PetscReal m1  = rho * h;
  const PetscReal m2  = rho * ftm;

  int                    ngps = 0;     // obere Schranke der Gauss-Punktschleife
  std::vector<PetscReal> weights;      // Integrationsgewichte


  // -------------------------------------------------------------------------
  //  Je nach Elementtyp (Quad oder Tria) werden die Integrationsgewichte
  //  zusammengesammelt.
  // -------------------------------------------------------------------------
  if (nnod == 4 || nnod == 9)
  {
    ngps = ngp * ngp;
    weights.resize(ngps);
    for (int n=0; n<ngps; n++)
      weights[n] = m_GaussPoints.getGaussWeight2D(ngp, n);
  }
  else if (nnod == 3 || nnod == 6)
  {
    ngps = ngp;
    weights.resize(ngps);
    for (int n=0; n<ngps; n++)
      weights[n] = 0.5 * m_GaussPoints.getGaussWeightTria(ngp, n);
  }


  // -------------------------------------------------------------------------
  //  Gauss-Punkt-Schleife
  // -------------------------------------------------------------------------
  for (int n=0;n<ngps;n++)
  {
    // -----------------------------------------------------------------------
    //  Aufstellen der Jacobi-Matrix
    // -----------------------------------------------------------------------
    setupJacobian2D(N, n);
    wdJ = weights[n] * detJac;

    // -----------------------------------------------------------------------
    //  Integration der Massenmatrix
    // -----------------------------------------------------------------------
    for (int z=0;z<nnod;z++)
    {
      for (int s=0;s<nnod;s++)
      {
        MM(z*ndofs  , s*ndofs  ) += N(0,z,n) * N(0,s,n) * m1 * wdJ;
        MM(z*ndofs+1, s*ndofs+1) += N(0,z,n) * N(0,s,n) * m2 * wdJ;
        MM(z*ndofs+2, s*ndofs+2) += N(0,z,n) * N(0,s,n) * m2 * wdJ;
      }
    }
  }
}



void cElementStructureMindlin::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  // -- wenn das Element nicht belastet ist, koennen wir hier aufhoeren
  if (m_ElementLoads.size() != 0)
  {

    const int nnod  = getNumberOfNodes();  // Anzahl der Knoten im Element
    const int ndofs = 3;                   // Anzahl der Freiheitsgrade je Knoten
    const short ngp = getNumberOfGaussPoints();
    cArray3d    N( getShapeFunctions() );

    int               ngps;        // obere Schranke der Gauss-Punkt-Schleife
    std::vector<PetscReal> weights;     // Integrationsgewichte

    double            wdJ;         // = Integrationsgewicht * det Jacobian

    cElementLoadStructure *ptrELS; // Zeiger auf eine Elementlast


    // -------------------------------------------------------------------------
    //  Je nach Elementtyp (Quad oder Tria) werden die Integrationsgewichte
    //  zusammengesammelt.
    // -------------------------------------------------------------------------
    if (nnod == 4 || nnod == 9)
    {
      ngps = ngp * ngp;
      weights.resize(ngps);
      for (int n=0; n<ngps; n++)
        weights[n] = m_GaussPoints.getGaussWeight2D(ngp, n);
    }
    else if (nnod == 3 || nnod == 6)
    {
      ngps = ngp;
      weights.resize(ngps);
      for (int n=0; n<ngps; n++)
        weights[n] = 0.5 * m_GaussPoints.getGaussWeightTria(ngp, n);
    }


    // -------------------------------------------------------------------------
    //  Schleife ueber alle Elementlasten
    // -------------------------------------------------------------------------
    std::multimap<short, cElementLoad *>::iterator it;
    for (it=m_ElementLoads.begin(); it!= m_ElementLoads.end(); it++)
    {
      ptrELS = dynamic_cast<cElementLoadStructure *>(it->second);

      if (ptrELS != NULL)
      {
        cElementVector Fz(nnod);
        for (int k=0; k<nnod; k++)
          Fz[k] = ptrELS->getForceComponent(m_Nodes[k]->getId(), 2);

        // -------------------------------------------------------------------------
        //  Gauss-Punkt-Schleife
        // -------------------------------------------------------------------------
        for (int n=0;n<ngp*ngp;n++)
        {
          // -----------------------------------------------------------------------
          //  Aufstellen der Jacobi-Matrix
          // -----------------------------------------------------------------------
          setupJacobian2D(N, n);

          // ---------------------------------------------------------------------
          //  Bestimmung des zum Integrationspunkt gehoerenden Gewichts
          // ---------------------------------------------------------------------
          wdJ = weights[n] * detJac;

          for (int z=0; z<nnod; z++)
            for (int s=0; s<nnod; s++)
              LV[z*ndofs] += N(0,z,n) * N(0,s,n) * Fz[s] * wdJ;
        }
      }
    }
    ptrELS = NULL;
  }
}
