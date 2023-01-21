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

#include "elementstructurebeam12.h"

#include <iomanip>

cElementStructureBeam12::cElementStructureBeam12(const eUseBeamTheory &Theory)
    : cElementStructureBeam3D(2, 6, 0) {
  m_UseTheory = Theory;
}

cElementStructureBeam12::cElementStructureBeam12(
    const cElementStructureBeam12 &other)
    : cElementStructureBeam3D(other) {
  m_UseTheory = other.getBeamTheory();
}

cElementStructureBeam12::~cElementStructureBeam12() {
  // empty
}

void cElementStructureBeam12::computeTransformationMatrix(
    cMatrix &T, const bool &transpose) const {
  // for(int i=0; i<12; i++)T(i,i)=1.0;

  // --- create transformationmatrix (12,12)
  //// --- we need a special transformation matrix in order
  ////     to get the right direction of the element
  cPoint a = *m_Nodes[0];
  cPoint b = *m_Nodes[1];

  cPoint new_a;
  cPoint new_b;
  for (int i = 0; i < 3; ++i) {
    if (a[i] < b[i]) {
      new_a[i] = a[i];
      new_b[i] = b[i];
    } else {
      new_a[i] = b[i];
      new_b[i] = a[i];
    }
  }

  // traslate notes
  cPoint p;
  for (int i = 0; i < 3; ++i) p[i] = new_b[i] - new_a[i];

  PetscReal l = sqrt(pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
  PetscReal r1 = sqrt(pow(p[0], 2) + pow(p[1], 2));

  double cos_beta = r1 / l;
  double sin_beta = p[2] / l;
  double cos_gamma, sin_gamma;

  if (r1 != 0.) {
    cos_gamma = p[0] / r1;
    sin_gamma = p[1] / r1;
  } else {
    cos_gamma = 0.;
    sin_gamma = 1.;
  }

  T(0, 0) = cos_beta * cos_gamma;
  T(0, 1) = cos_beta * sin_gamma;
  T(0, 2) = sin_beta;
  T(1, 0) = -sin_gamma;
  T(1, 1) = cos_gamma;
  T(1, 2) = 0.0;
  T(2, 0) = -sin_beta * cos_gamma;
  T(2, 1) = -sin_beta * sin_gamma;
  T(2, 2) = cos_beta;

  T(3, 3) = cos_beta * cos_gamma;
  T(3, 4) = cos_beta * sin_gamma;
  T(3, 5) = sin_beta;
  T(4, 3) = -sin_gamma;
  T(4, 4) = cos_gamma;
  T(4, 5) = 0.0;
  T(5, 3) = -sin_beta * cos_gamma;
  T(5, 4) = -sin_beta * sin_gamma;
  T(5, 5) = cos_beta;

  T(6, 6) = cos_beta * cos_gamma;
  T(6, 7) = cos_beta * sin_gamma;
  T(6, 8) = sin_beta;
  T(7, 6) = -sin_gamma;
  T(7, 7) = cos_gamma;
  T(7, 8) = 0.0;
  T(8, 6) = -sin_beta * cos_gamma;
  T(8, 7) = -sin_beta * sin_gamma;
  T(8, 8) = cos_beta;

  T(9, 9) = cos_beta * cos_gamma;
  T(9, 10) = cos_beta * sin_gamma;
  T(9, 11) = sin_beta;
  T(10, 9) = -sin_gamma;
  T(10, 10) = cos_gamma;
  T(10, 11) = 0.0;
  T(11, 9) = -sin_beta * cos_gamma;
  T(11, 10) = -sin_beta * sin_gamma;
  T(11, 11) = cos_beta;

  if (transpose == true) {
    infam::transpose(T, T);
  }
}

std::vector<eKnownDofs> cElementStructureBeam12::getDofs(void) const {
  std::vector<eKnownDofs> res(6);

  res[0] = disp_x1;
  res[1] = disp_x2;
  res[2] = disp_x3;
  res[3] = disp_w1;
  res[4] = disp_w2;
  res[5] = disp_w3;
  return res;
}

void cElementStructureBeam12::assembleStiffnessMatrix(cElementMatrix &KM,
                                                      Vec *x = NULL,
                                                      Vec *dx = NULL) {
  cElementMatrix K(12, 12);
  PetscScalar E;

  const PetscReal A = m_Material->getA();
  const PetscReal nu = m_Material->getNu();
  const PetscReal Ix = m_Material->getIx();
  const PetscReal Iy = m_Material->getIy();
  const PetscReal Iz = m_Material->getIz();
  const PetscReal L = getLength();

  E = m_Material->getEOmega();

  const PetscReal PSI = getPsi();
  PetscReal etaY;
  PetscReal etaZ;

  if (PSI != 1.0) {
    etaY = PSI * Iy;
    etaZ = PSI * Iz;
  } else {
    etaY = 0.0;
    etaZ = 0.0;
  }
  const PetscScalar G = E / (2. * (1. + nu));

  const PetscScalar k1 = E * A / L;
  const PetscScalar k2 = 12. * E * Iz / (L * L * L) * 1. / (1. + etaY);
  const PetscScalar k3 = 6. * E * Iz / (L * L) * 1. / (1. + etaY);
  const PetscScalar k4 = E * Iz / L * (4. + etaY) / (1. + etaY);
  const PetscScalar k5 = E * Iz / L * (2. - etaY) / (1. + etaY);
  const PetscScalar k6 = 12. * E * Iy / (L * L * L) * 1. / (1. + etaZ);
  const PetscScalar k7 = 6. * E * Iy / (L * L) * 1. / (1. + etaZ);
  const PetscScalar k8 = E * Iy / L * (4. + etaZ) / (1. + etaZ);
  const PetscScalar k9 = E * Iy / L * (2. - etaZ) / (1. + etaZ);
  const PetscScalar k10 = G * Ix / L;

  K(0, 0) = K(6, 6) = k1;
  K(0, 6) = K(6, 0) = -k1;

  K(3, 3) = K(9, 9) = k10;
  K(3, 9) = K(9, 3) = -k10;

  K(2, 2) = K(8, 8) = k6;
  K(8, 2) = K(2, 8) = -k6;
  K(8, 4) = K(4, 8) = k7;
  K(8, 10) = K(10, 8) = k7;
  K(2, 4) = K(4, 2) = -k7;
  K(2, 10) = K(10, 2) = -k7;
  K(4, 4) = K(10, 10) = k8;
  K(4, 10) = K(10, 4) = k9;

  K(1, 1) = K(7, 7) = k2;
  K(1, 7) = K(7, 1) = -k2;
  K(1, 5) = K(5, 1) = k3;
  K(1, 11) = K(11, 1) = k3;
  K(5, 7) = K(7, 5) = -k3;
  K(7, 11) = K(11, 7) = -k3;
  K(5, 5) = K(11, 11) = k4;
  K(5, 11) = K(11, 5) = k5;

  // --- compute KM = T^t * K * T
  cMatrix T(12, 12);
  computeTransformationMatrix(T);
  BT_C_B_wdJ(KM, T, K, 1.);
}

void cElementStructureBeam12::assembleMassMatrix(cElementMatrix &MM) {
  cElementMatrix M(12, 12);

  const PetscReal A = m_Material->getA();
  const PetscReal Ix = m_Material->getIx();
  // const PetscReal   Iy    = m_Material->getIy();
  // const PetscReal   Iz    = m_Material->getIz();
  const PetscReal rho = m_Material->getRho();
  const PetscReal L = getLength();
  // const PetscReal   PSI   = getPsi();

  const PetscReal m = rho * A * L;  // Mass
  const PetscReal m1 = one_over_three * m;
  const PetscReal m2 = one_over_six * m;
  const PetscReal m3 = 13. / 35. * m;
  const PetscReal m4 = 11. / 210. * m * L;
  const PetscReal m5 = 9. / 70. * m;
  const PetscReal m6 = 13. / 420. * m * L;
  const PetscReal m7 = 1. / 105. * m * L * L;
  const PetscReal m8 = 1. / 140. * m * L * L;
  const PetscReal m9 = one_over_three * rho * Ix * L;
  const PetscReal m10 = one_over_six * rho * Ix * L;

  M(0, 0) = M(6, 6) = m1;
  M(0, 6) = M(6, 0) = m2;

  M(3, 3) = M(9, 9) = m9;
  M(3, 9) = M(9, 3) = m10;

  M(1, 1) = M(2, 2) = m3;
  M(7, 7) = M(8, 8) = m3;
  M(4, 4) = M(5, 5) = m7;
  M(10, 10) = M(11, 11) = m7;

  M(1, 5) = M(5, 1) = m4;
  M(8, 10) = M(10, 8) = m4;
  M(2, 4) = M(4, 2) = -m4;
  M(7, 11) = M(11, 7) = -m4;

  M(1, 7) = M(7, 1) = m5;
  M(2, 8) = M(8, 2) = m5;
  M(4, 10) = M(10, 4) = -m8;
  M(5, 11) = M(11, 5) = -m8;

  M(1, 11) = M(11, 1) = -m6;
  M(2, 10) = M(10, 2) = m6;
  M(4, 8) = M(8, 4) = -m6;
  M(5, 7) = M(7, 5) = m6;

  // --- compute MM = T^t * M * T
  cMatrix T(12, 12);
  computeTransformationMatrix(T);
  BT_C_B_wdJ(MM, T, M, 1.);
}

void cElementStructureBeam12::assembleLoadVector(cElementVector &LV,
                                                 cElementMatrix &KM,
                                                 Vec *x = NULL,
                                                 Vec *dx = NULL) {
  // --- leave if no load is applied
  //     to this element
  if (m_ElementLoads.size() != 0) {
    // --- Px and Py are constant elementloads whose x and y direction
    //     correspond to the GLOBAL coordinate system
    PetscScalar Px = 0.;
    PetscScalar Py = 0.;
    PetscScalar Pz = 0.;
    PetscScalar Mx = 0.;
    const PetscReal L = getLength();

    // --- first the transformation matrix T is used
    //     to compute the local forces of this element
    cMatrix T(12, 12);
    computeTransformationMatrix(T);
    // --- load vector in LOCAL coordinate system
    cElementVector LVL(12);

    std::multimap<short, cElementLoad *>::iterator it;
    cElementLoadStructure *ptrELS = NULL;

    for (it = m_ElementLoads.begin(); it != m_ElementLoads.end(); it++) {
      ptrELS = dynamic_cast<cElementLoadStructure *>(it->second);

      if (ptrELS != NULL) {
        // --- we assume the load to be constant
        Px = ptrELS->getForceComponent(m_Nodes[0]->getId(), 0);
        Py = ptrELS->getForceComponent(m_Nodes[0]->getId(), 1);
        Pz = ptrELS->getForceComponent(m_Nodes[0]->getId(), 2);

        // --- we need the force components in local
        //     coordinate directions
        // std::cout<<Px<<" "<<Py<<" "<<" "<<Pz<<std::endl;
        const PetscScalar PLx = T(0, 0) * Px + T(0, 1) * Py + T(0, 2) * Pz;
        const PetscScalar PLy = T(1, 0) * Px + T(1, 1) * Py + T(1, 2) * Pz;
        const PetscScalar PLz = T(2, 0) * Px + T(2, 1) * Py + T(2, 2) * Pz;
        // std::cout<<PLx<<" "<<PLy<<" "<<" "<<PLz<<std::endl;

        LVL[0] += 1. / 2. * PLx * L;
        LVL[1] += 1. / 2. * PLy * L;
        LVL[2] += 1. / 2. * PLz * L;
        LVL[3] += 1. / 2. * Mx * L;
        LVL[4] += -1. / 12. * PLz * L * L;
        LVL[5] += 1. / 12. * PLy * L * L;
        LVL[6] += 1. / 2. * PLx * L;
        LVL[7] += 1. / 2. * PLy * L;
        LVL[8] += 1. / 2. * PLz * L;
        LVL[9] += 1. / 2. * Mx * L;
        LVL[10] += 1. / 12. * PLz * L * L;
        LVL[11] += -1. / 12. * PLy * L * L;
      }
      ptrELS = NULL;

      // --- setup transformationmatrix
      //     and compute LV += T^t * LVL

      computeTransformationMatrix(T, true);  // now we need T^t
      infam::mult(T, LVL, LV);
    }
  }
  assembleInternalForces(x, KM, LV);
}

std::ostream &cElementStructureBeam12::write(std::ostream &os) const {
  if (m_UseTheory == Bernoulli)
    os << "beam3D element (Euler-Bernoulli3D) 12dofs" << std::endl;
  else
    os << "beam3D element (Timoshenko3D) 12dofs" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k = 0; k < getNumberOfNodes(); k++) os << *(m_Nodes[k]);
  return os;
}

std::ostream &cElementStructureBeam12::writeXml(std::ostream &os) const {
  if (m_UseTheory == Bernoulli)
    os << "<BeamBernoulli12>";
  else
    os << "<BeamTimoshenko12>";
  os << "<Id>" << getId() << "</Id>";
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  if (m_UseTheory == Bernoulli)
    os << "</BeamBernoulli12>";
  else
    os << "</BeamTimoshenko12>";
  return os;
}
