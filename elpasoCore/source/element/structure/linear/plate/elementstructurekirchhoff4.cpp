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

#include "elementstructurekirchhoff4.h"


cElementStructureKirchhoff4::cElementStructureKirchhoff4() :
  cElementStructureKirchhoff(4, 4)
{
  if (howMany() == 1)
    initializeShapeFunctions();
}


cElementStructureKirchhoff4::cElementStructureKirchhoff4(const cElementStructureKirchhoff4 &other) :
  cElementStructureKirchhoff(other)
{
  // empty
}


cElementStructureKirchhoff4::~cElementStructureKirchhoff4()
{
  // empty
}


std::ostream& cElementStructureKirchhoff4::write(std::ostream &os) const
{
  os << "Kirchhoff plateelement no. " << getId() << std::endl;
  os << "  number of nodes : " << getNumberOfNodes() << std::endl;
  os << "  degrees of freedom per node : " << getNumberOfDofsPerNode() << std::endl;

  for (int i=0;i<getNumberOfNodes();i++)
    os << *m_Nodes[i] << std::endl;

  os << *m_Material;

  return os;
}


// ---------------------------------------------------------------------------
//  initialize static class members
// ---------------------------------------------------------------------------
cArray3d cElementStructureKirchhoff4::N;


void cElementStructureKirchhoff4::initializeShapeFunctions(void)
{
  const int nnod = getNumberOfNodes();
  const int ngp  = getNumberOfGaussPoints();

  // -------------------------------------------------------------------------
  //  allocate memory
  // -------------------------------------------------------------------------
  N.initialize(3,4,ngp*ngp);


  // -------------------------------------------------------------------------
  //  initialize test functions
  // -------------------------------------------------------------------------
  for (int n=0; n<ngp*ngp; n++)
  {
    cPoint gp( m_GaussPoints.getGaussPoint2D(ngp, n) );

    for (int k=0; k<nnod; k++)
    {
      N(0,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_fun, gp);
      N(1,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_xi , gp);
      N(2,k,n) = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k, N_eta, gp);
    }
  }
}



void cElementStructureKirchhoff4::assembleStiffnessMatrix(cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  PetscReal wdJ;

  cPoint    gp;
  PetscReal funkt11,funkt21,funkt31,funkt41;
  PetscReal funkt12,funkt22,funkt32,funkt42;
  PetscReal ablksi2[16], ableta2[16], ablksi1eta1[16];

  rotateElement();


  // --- elasticitymatrix
  cElementMatrix Cb(3,3);
  m_Material->setupCb(Cb, this);
  
  // --- Gauss point loop
  const int ngp = getNumberOfGaussPoints();
  for (int n=0;n<ngp*ngp;n++)
  {
    // --- coordinates of current Gauss point
    gp = m_GaussPoints.getGaussPoint2D(ngp, n);

    // --- compute Jacobian
    setupJacobian2D(N, n);

    const PetscReal lx = Jac(0,0); // side lengths of element
    const PetscReal ly = Jac(1,1);

    // --- add contribution of current Gauss point
    wdJ = m_GaussPoints.getGaussWeight2D(ngp, n);

    const PetscScalar ed11 = detJac/lx/lx * Cb(0,0) /lx/lx   ;
    const PetscScalar ed12 = detJac/lx/lx * Cb(0,1) /ly/ly   ;
    const PetscScalar ed21 = detJac/ly/ly * Cb(1,0) /lx/lx   ;
    const PetscScalar ed22 = detJac/ly/ly * Cb(1,1) /ly/ly   ;
    const PetscScalar ed33 = detJac/lx/ly * Cb(2,2) /lx/ly   ;

    const PetscReal F1 = gp[0];
    const PetscReal F2 = gp[1];
    funkt11 =  1.5*F1;
    funkt21 = (-0.5+1.5*F1) * lx;
    funkt31 = -1.5 * F1;
    funkt41 = ( 0.5+1.5*F1) * lx;
    funkt12 = (2.0-3.0*F2+F2*F2*F2) / 4.0;
    funkt22 = ( 1.0-F2-F2*F2+F2*F2*F2) * ly / 4.0;
    funkt32 = (2.0+3.0*F2-F2*F2*F2) / 4.0;
    funkt42 = (-1.0-F2+F2*F2+F2*F2*F2) * ly / 4.0;
    ablksi2[ 0] = funkt11 * funkt12;
    ablksi2[ 1] = funkt21 * funkt12;
    ablksi2[ 2] = funkt11 * funkt22;
    ablksi2[ 3] = funkt21 * funkt22;
    ablksi2[ 4] = funkt31 * funkt12;
    ablksi2[ 5] = funkt41 * funkt12;
    ablksi2[ 6] = funkt31 * funkt22;
    ablksi2[ 7] = funkt41 * funkt22;
    ablksi2[ 8] = funkt31 * funkt32;
    ablksi2[ 9] = funkt41 * funkt32;
    ablksi2[10] = funkt31 * funkt42;
    ablksi2[11] = funkt41 * funkt42;
    ablksi2[12] = funkt11 * funkt32;
    ablksi2[13] = funkt21 * funkt32;
    ablksi2[14] = funkt11 * funkt42;
    ablksi2[15] = funkt21 * funkt42;
    funkt11 = (2.0-3.0*F1+F1*F1*F1) / 4.0;
    funkt21 = ( 1.0-F1-F1*F1+F1*F1*F1) * lx / 4.0;
    funkt31 = (2.0+3.0*F1-F1*F1*F1) / 4.0;
    funkt41 = (-1.0-F1+F1*F1+F1*F1*F1) * lx / 4.0;
    funkt12 =  1.5 * F2     ;
    funkt22 = (-0.5+1.5*F2) * ly  ;
    funkt32 = -1.5 * F2     ;
    funkt42 = ( 0.5+1.5*F2) * ly  ;
    ableta2[ 0] = funkt11 * funkt12;
    ableta2[ 1] = funkt21 * funkt12;
    ableta2[ 2] = funkt11 * funkt22;
    ableta2[ 3] = funkt21 * funkt22;
    ableta2[ 4] = funkt31 * funkt12;
    ableta2[ 5] = funkt41 * funkt12;
    ableta2[ 6] = funkt31 * funkt22;
    ableta2[ 7] = funkt41 * funkt22;
    ableta2[ 8] = funkt31 * funkt32;
    ableta2[ 9] = funkt41 * funkt32;
    ableta2[10] = funkt31 * funkt42;
    ableta2[11] = funkt41 * funkt42;
    ableta2[12] = funkt11 * funkt32;
    ableta2[13] = funkt21 * funkt32;
    ableta2[14] = funkt11 * funkt42;
    ableta2[15] = funkt21 * funkt42;
    funkt11 = -0.75 + 0.75 * F1 * F1;
    funkt21 = (-1.0-2.0*F1+3.0*F1*F1) * lx / 4.0;
    funkt31 =  0.75 - 0.75 * F1 * F1;
    funkt41 = (-1.0+2.0*F1+3.0*F1*F1) * lx / 4.0;
    funkt12 = -0.75 + 0.75 * F2 * F2;
    funkt22 = (-1.0-2.0*F2+3.0*F2*F2) * ly / 4.0;
    funkt32 =  0.75 - 0.75 * F2 * F2;
    funkt42 = (-1.0+2.0*F2+3.0*F2*F2) * ly / 4.0;
    ablksi1eta1[ 0] = funkt11 * funkt12;
    ablksi1eta1[ 1] = funkt21 * funkt12;
    ablksi1eta1[ 2] = funkt11 * funkt22;
    ablksi1eta1[ 3] = funkt21 * funkt22;
    ablksi1eta1[ 4] = funkt31 * funkt12;
    ablksi1eta1[ 5] = funkt41 * funkt12;
    ablksi1eta1[ 6] = funkt31 * funkt22;
    ablksi1eta1[ 7] = funkt41 * funkt22;
    ablksi1eta1[ 8] = funkt31 * funkt32;
    ablksi1eta1[ 9] = funkt41 * funkt32;
    ablksi1eta1[10] = funkt31 * funkt42;
    ablksi1eta1[11] = funkt41 * funkt42;
    ablksi1eta1[12] = funkt11 * funkt32;
    ablksi1eta1[13] = funkt21 * funkt32;
    ablksi1eta1[14] = funkt11 * funkt42;
    ablksi1eta1[15] = funkt21 * funkt42;

    // ke-Matrix berechnen
    for(int kz=0; kz<16; kz++)
    for(int ks=0; ks<16; ks++)
      KM(kz,ks) += wdJ *
      ( ablksi2[kz] * ed11 * ablksi2[ks]
      + ablksi2[kz] * ed12 * ableta2[ks]
      + ableta2[kz] * ed21 * ablksi2[ks]
      + ableta2[kz] * ed22 * ableta2[ks]
      + ablksi1eta1[kz] * ed33 * ablksi1eta1[ks] * 4.0);
  }
}


void cElementStructureKirchhoff4::assembleMassMatrix(cElementMatrix &MM)
{
  PetscReal wdJ;
  cPoint    gp;
  PetscReal funkt11,funkt21,funkt31,funkt41;
  PetscReal funkt12,funkt22,funkt32,funkt42;
  PetscReal bikub[16];
  PetscReal F1, F2;
  PetscReal lx, ly;


  // --- used materialparameters
  const PetscReal m1  = m_Material->getRho() * m_Material->getT(this);

  rotateElement();

  // --- Gauss point loop
  const int ngp = getNumberOfGaussPoints();
  const int nnod = getNumberOfNodes();
  const int ndofs = getNumberOfDofsPerNode();

  for (int n=0;n<ngp*ngp;n++)
  {
    // --- compute Jacobian
    setupJacobian2D(N, n);

    gp  = m_GaussPoints.getGaussPoint2D(ngp, n);
    wdJ = m_GaussPoints.getGaussWeight2D(ngp, n) * detJac;

    lx = Jac(0,0);
    ly = Jac(1,1);

    F1 = gp[0];
    F2 = gp[1];
    funkt11 = 0.25 * (2.0 - 3.0*F1 + F1*F1*F1); // Herm[1]
    funkt12 = 0.25 * (2.0 - 3.0*F2 + F2*F2*F2);
    funkt31 = 0.25 * (2.0 + 3.0*F1 - F1*F1*F1); // Herm[2]
    funkt32 = 0.25 * (2.0 + 3.0*F2 - F2*F2*F2);
    funkt21 = ( 1.0 - F1 - F1*F1 + F1*F1*F1) * 0.25 * lx;  // Herm[3]
    funkt22 = ( 1.0 - F2 - F2*F2 + F2*F2*F2) * 0.25 * ly;
    funkt41 = (-1.0 - F1 + F1*F1 + F1*F1*F1) * 0.25 * lx;  // Herm[4]
    funkt42 = (-1.0 - F2 + F2*F2 + F2*F2*F2) * 0.25 * ly;
    bikub[ 0] = funkt11 * funkt12;
    bikub[ 1] = funkt21 * funkt12;
    bikub[ 2] = funkt11 * funkt22;
    bikub[ 3] = funkt21 * funkt22;
    bikub[ 4] = funkt31 * funkt12;
    bikub[ 5] = funkt41 * funkt12;
    bikub[ 6] = funkt31 * funkt22;
    bikub[ 7] = funkt41 * funkt22;
    bikub[ 8] = funkt31 * funkt32;
    bikub[ 9] = funkt41 * funkt32;
    bikub[10] = funkt31 * funkt42;
    bikub[11] = funkt41 * funkt42;
    bikub[12] = funkt11 * funkt32;
    bikub[13] = funkt21 * funkt32;
    bikub[14] = funkt11 * funkt42;
    bikub[15] = funkt21 * funkt42;

    // --- compute massmatrix
    for (int z=0;z<nnod*ndofs; z++)
    for (int s=0;s<nnod*ndofs; s++)
      MM(z,s) += bikub[z] * bikub[s] * wdJ * m1;
  }
}


void cElementStructureKirchhoff4::assembleLoadVector(cElementVector &LV, cElementMatrix &KM, Vec *x = NULL, Vec *dx = NULL)
{
  PetscReal F1,F2,fak;
  PetscReal funkt11,funkt21,funkt31,funkt41;
  PetscReal funkt12,funkt22,funkt32,funkt42;
  PetscReal bikub[16];
  PetscReal lx, ly;
  cPoint gp;
  cElementLoadStructure *last = NULL;
  std::multimap<short, cElementLoad *>::iterator it;
  const int ngp = getNumberOfGaussPoints();

  cElementVector Fz(4);


  for (it=m_ElementLoads.begin(); it!=m_ElementLoads.end(); it++)
  {
    last = dynamic_cast<cElementLoadStructure *>(it->second);

    if (last != NULL)
    {
      for (int k=0; k<4; k++)
        Fz[k] = last->getForceComponent(m_Nodes[k]->getId(), 2);
      last = NULL;

      for(int n=0; n<ngp*ngp; n++)
      {
        gp = m_GaussPoints.getGaussPoint2D(ngp, n);

        // --- compute Jacobian
        setupJacobian2D(N, n);
        lx = Jac(0,0);
        ly = Jac(1,1);


        F1 = gp[0];
        F2 = gp[1];
        funkt11 = 0.25 * (2.0-3.0*F1+F1*F1*F1);
        funkt21 = ( 1.0-F1-F1*F1+F1*F1*F1) * 0.25 * lx;
        funkt31 = 0.25 * (2.0+3.0*F1-F1*F1*F1);
        funkt41 = (-1.0-F1+F1*F1+F1*F1*F1) * 0.25 * lx;
        funkt12 = 0.25 * (2.0-3.0*F2+F2*F2*F2);
        funkt22 = ( 1.0-F2-F2*F2+F2*F2*F2) * 0.25 * ly;
        funkt32 = 0.25 * (2.0+3.0*F2-F2*F2*F2);
        funkt42 = (-1.0-F2+F2*F2+F2*F2*F2) * 0.25 * ly;
        bikub[ 0] = funkt11 * funkt12;
        bikub[ 1] = funkt21 * funkt12;
        bikub[ 2] = funkt11 * funkt22;
        bikub[ 3] = funkt21 * funkt22;
        bikub[ 4] = funkt31 * funkt12;
        bikub[ 5] = funkt41 * funkt12;
        bikub[ 6] = funkt31 * funkt22;
        bikub[ 7] = funkt41 * funkt22;
        bikub[ 8] = funkt31 * funkt32;
        bikub[ 9] = funkt41 * funkt32;
        bikub[10] = funkt31 * funkt42;
        bikub[11] = funkt41 * funkt42;
        bikub[12] = funkt11 * funkt32;
        bikub[13] = funkt21 * funkt32;
        bikub[14] = funkt11 * funkt42;
        bikub[15] = funkt21 * funkt42;

        // --- compute elementload
        fak = m_GaussPoints.getGaussWeight2D(ngp, n) * lx * ly;

        for(int k=0; k<16; k++)
        for(int j=0; j<4;  j++)
          LV[k] += bikub[k] * N(0,j,n) * fak * Fz[j];
      }
    }
  }
}


// ----------------------------------------------------------------------------
// --- check element for proper alignment
// ----------------------------------------------------------------------------
void cElementStructureKirchhoff4::rotateElement()
{
  PetscReal lx, ly;

  int counting = 0;

  for(;;)
  {
    lx = (*m_Nodes[1])[0] - (*m_Nodes[0])[0];
    ly = (*m_Nodes[2])[1] - (*m_Nodes[0])[1];

    if ( (lx < cstGeomEps) || (ly < cstGeomEps) )
    {
      counting++;
      std::rotate(m_Nodes.begin(), m_Nodes.begin() + 1, m_Nodes.end());

      if (counting > 4)
        throw cException("inappropriate Kirchhoff plateelement!", __FILE__, __LINE__);
    }
    else
      return;
  }
}




std::ostream& cElementStructureKirchhoff4::writeXml(std::ostream &os) const
{
  os << "<Kirch4>";
  os << "<Id>" << getId() << "</Id>";
  for (int k=0; k<getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</Kirch4>";
  return os;
}
