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

#include "elementstructuremindlindsg4.h"

cElementStructureMindlinDSG4::cElementStructureMindlinDSG4(void)
    : cElementStructureMindlin(4, 2) {
  // empty
}

cElementStructureMindlinDSG4::cElementStructureMindlinDSG4(
    const cElementStructureMindlinDSG4 &other)
    : cElementStructureMindlin(other) {
  // empty
}

cElementStructureMindlinDSG4::~cElementStructureMindlinDSG4() {
  // empty
}

void cElementStructureMindlinDSG4::assembleStiffnessMatrix(cElementMatrix &KM,
                                                           Vec *x = NULL,
                                                           Vec *dx = NULL) {
  PetscReal wdJ;  // Gewicht * detJ
  const int nnod = getNumberOfNodes();
  const short ngp = getNumberOfGaussPoints();
  const int ndofs = 3;
  cMatrix Bb(3, nnod * ndofs);  // bending part
  cMatrix Bs(2, nnod * ndofs);  // shear part

  // -------------------------------------------------------------------------
  //   Stabilization according to Lyly and Stenberg
  // -------------------------------------------------------------------------
  PetscReal hk = 0.0;            // typical element size
  const PetscReal alpha = 0.12;  // weighting factor
  int sequenz[] = {0, 1, 2, 3, 0};
  const PetscReal thick = m_Material->getT(this);

  // -------------------------------------------------------------------------
  //   Fuer hk sind alle Werte moeglich, die die Groesse des Elements
  //   angemessen repraesentieren, wie z.B. die maximale Seitenlaenge
  //   (wird hier verwendet) oder die Wurzel der Elementflaeche
  // -------------------------------------------------------------------------
  for (int i = 1; i < 4; i++)
    hk = std::max(hk,
                  (*m_Nodes[sequenz[i]]).distance((*m_Nodes[sequenz[i - 1]])));

  const PetscReal tau = 1.0 / (1.0 + alpha * hk * hk / (thick * thick));

  cElementMatrix Cs(2, 2);
  cElementMatrix Cb(3, 3);

  m_Material->setupCs(Cs, this);
  m_Material->setupCb(Cb, this);

  infam::scale(Cs, tau);

  // Check orientation (global or local) and rotate material matrix if necessary
  if (m_Orientation ==
      "local")  // Material matrix is given in global system -> Transformation
                // needed for local coordinate system
  {
    // Local x direction in global coordinate system
    PetscReal lx =
        0.5 * (m_Nodes[1]->getComponent(0) + m_Nodes[2]->getComponent(0)) -
        0.5 * (m_Nodes[0]->getComponent(0) + m_Nodes[3]->getComponent(0));
    PetscReal ly =
        0.5 * (m_Nodes[1]->getComponent(1) + m_Nodes[2]->getComponent(1)) -
        0.5 * (m_Nodes[0]->getComponent(1) + m_Nodes[3]->getComponent(1));
    PetscReal vectorLengthl = std::sqrt(lx * lx + ly * ly);
    lx = lx / vectorLengthl;
    ly = ly / vectorLengthl;
    // Local y direction in global coordinate system is calculated
    PetscReal rx = -ly;
    PetscReal ry = lx;

    cMatrix Gb(3,
               3);  // Rotation matrix G for bending part (Cb material matrix)
    Gb(0, 0) = lx * lx;
    Gb(0, 1) = ly * ly;
    Gb(0, 2) = lx * ly;
    Gb(1, 0) = rx * rx;
    Gb(1, 1) = ry * ry;
    Gb(1, 2) = rx * ry;
    Gb(2, 0) = 2 * lx * rx;
    Gb(2, 1) = 2 * ly * ry;
    Gb(2, 2) = lx * ry + ly * rx;

    cElementMatrix CbTemp(3, 3);
    infam::BT_C_B_wdJ(CbTemp, Gb, Cb, 1.0);  // Rotation of material matrix
    Cb = CbTemp;
  }
  if (m_Orientation ==
      "user-def")  // Material matrix is given in global system ->
                   // Transformation needed for user-defined coordinate system
  {
    PetscReal lx = m_OrientationVector[0];
    PetscReal ly = m_OrientationVector[1];

    PetscReal rx = m_OrientationVector[3];
    PetscReal ry = m_OrientationVector[4];

    cMatrix Gb(3,
               3);  // Rotation matrix G for bending part (Cb material matrix)
    Gb(0, 0) = lx * lx;
    Gb(0, 1) = ly * ly;
    Gb(0, 2) = lx * ly;
    Gb(1, 0) = rx * rx;
    Gb(1, 1) = ry * ry;
    Gb(1, 2) = rx * ry;
    Gb(2, 0) = 2 * lx * rx;
    Gb(2, 1) = 2 * ly * ry;
    Gb(2, 2) = lx * ry + ly * rx;

    cElementMatrix CbTemp(3, 3);
    infam::BT_C_B_wdJ(CbTemp, Gb, Cb, 1.0);  // Rotation of material matrix
    Cb = CbTemp;
  }

  PetscReal dNdx, dNdy;

  // -------------------------------------------------------------------------
  //  Gauss point loop
  // -------------------------------------------------------------------------
  for (int n = 0; n < ngp * ngp; n++) {
    // -----------------------------------------------------------------------
    //  setup Jacobian
    // -----------------------------------------------------------------------
    setupJacobian2D(N, n);

    const PetscReal a = Jac(0, 0);
    const PetscReal b = Jac(0, 1);
    const PetscReal d = Jac(1, 0);
    const PetscReal c = Jac(1, 1);

    wdJ = m_GaussPoints.getGaussWeight2D(ngp, n) * detJac;

    // -----------------------------------------------------------------------
    //   setup Bb - bending
    // -----------------------------------------------------------------------
    for (int k = 0; k < nnod; k++) {
      dNdx = (c * N(1, k, n) - b * N(2, k, n)) / detJac;
      dNdy = (-d * N(1, k, n) + a * N(2, k, n)) / detJac;

      Bb(0, k * ndofs + 2) = dNdx;
      Bb(1, k * ndofs + 1) = -dNdy;
      Bb(2, k * ndofs + 2) = dNdy;
      Bb(2, k * ndofs + 1) = -dNdx;
    }

    // -----------------------------------------------------------------------
    //   setup Bs - shear
    // -----------------------------------------------------------------------
    Bs.setValue(0.0);
    Bs(0, 0) = -c * N(1, 1, n) + b * N(2, 3, n);  // node 1
    Bs(0, 1) = -b * c * N(1, 1, n) + b * c * N(2, 3, n);
    Bs(0, 2) = a * c * N(1, 1, n) - b * d * N(2, 3, n);

    Bs(1, 0) = d * N(1, 1, n) - a * N(2, 3, n);
    Bs(1, 1) = -a * c * N(2, 3, n) + b * d * N(1, 1, n);
    Bs(1, 2) = a * d * N(2, 3, n) - a * d * N(1, 1, n);

    Bs(0, 3) = c * N(1, 1, n) + b * N(2, 2, n);  // node 2
    Bs(0, 4) = -b * c * N(1, 1, n) + b * c * N(2, 2, n);
    Bs(0, 5) = a * c * N(1, 1, n) - b * d * N(2, 2, n);

    Bs(1, 3) = -a * N(2, 2, n) - d * N(1, 1, n);
    Bs(1, 4) = -a * c * N(2, 2, n) + b * d * N(1, 1, n);
    Bs(1, 5) = a * d * N(2, 2, n) - a * d * N(1, 1, n);

    Bs(0, 6) = c * N(1, 2, n) - b * N(2, 2, n);  // node 3
    Bs(0, 7) = -b * c * N(1, 2, n) + b * c * N(2, 2, n);
    Bs(0, 8) = a * c * N(1, 2, n) - b * d * N(2, 2, n);

    Bs(1, 6) = a * N(2, 2, n) - d * N(1, 2, n);
    Bs(1, 7) = -a * c * N(2, 2, n) + b * d * N(1, 2, n);
    Bs(1, 8) = a * d * N(2, 2, n) - a * d * N(1, 2, n);

    Bs(0, 9) = -c * N(1, 2, n) - b * N(2, 3, n);  // node 4
    Bs(0, 10) = -b * c * N(1, 2, n) + b * c * N(2, 3, n);
    Bs(0, 11) = a * c * N(1, 2, n) - b * d * N(2, 3, n);

    Bs(1, 9) = a * N(2, 3, n) + d * N(1, 2, n);
    Bs(1, 10) = -a * c * N(2, 3, n) + b * d * N(1, 2, n);
    Bs(1, 11) = a * d * N(2, 3, n) - a * d * N(1, 2, n);

    infam::scale(Bs, 1. / detJac);

    // -----------------------------------------------------------------------
    //   integrate
    // -----------------------------------------------------------------------
    infam::BT_C_B_wdJ(KM, Bb, Cb, wdJ);
    infam::BT_C_B_wdJ(KM, Bs, Cs, wdJ);
  }
}

std::ostream &cElementStructureMindlinDSG4::write(std::ostream &os) const {
  os << "DSG4 element" << std::endl;
  os << "Id = " << getId() << std::endl;
  os << *m_Material << std::endl;

  for (int k = 0; k < getNumberOfNodes(); k++) os << *(m_Nodes[k]) << std::endl;

  return os;
}

std::ostream &cElementStructureMindlinDSG4::writeXml(std::ostream &os) const {
  os << "<DSG4>";
  os << "<Id>" << getId() << "</Id>";
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</DSG4>";
  return os;
}
