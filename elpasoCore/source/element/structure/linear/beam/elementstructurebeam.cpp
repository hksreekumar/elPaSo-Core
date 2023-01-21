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

#include "elementstructurebeam.h"

cElementStructureBeam::cElementStructureBeam(const eUseBeamTheory &Theory)
    : cElementStructureLinear(2, 3, 0) {
  m_UseTheory = Theory;
}

cElementStructureBeam::cElementStructureBeam(const cElementStructureBeam &other)
    : cElementStructureLinear(other) {
  m_UseTheory = other.getBeamTheory();
}

cElementStructureBeam::~cElementStructureBeam() {
  // empty
}

cVector cElementStructureBeam::getGlobalNormalVector(void) const {
  const PetscReal L = getLength();
  PetscReal nx = ((*m_Nodes[1])[0] - (*m_Nodes[0])[0]) / L;
  PetscReal ny = ((*m_Nodes[1])[1] - (*m_Nodes[0])[1]) / L;

  // --- we need a special transformation matrix in order
  //     to get the right direction of the element
  if (nx < 0.) {
    nx *= -1.;
    ny *= -1.;
  } else if ((nx < cstGeomEps) && (ny < 0.)) {
    nx *= -1.;
    ny *= -1.;
  }

  // --- local normal vector n^T = [ 0, 1, 0 ];
  //     globalNormal = T^t * n
  cVector globalNormal(3);

  globalNormal[0] = -ny;  // T^t(0,1)
  globalNormal[1] = nx;   // T^t(1,1)
  globalNormal[2] = 0.;

  infam::scale(globalNormal, 1. / globalNormal.abs());

  return globalNormal;
}

void cElementStructureBeam::computeTransformationMatrix(
    cMatrix &T, const bool &transpose) const {
  const PetscReal L = getLength();
  const PetscReal dx = (*m_Nodes[1])[0] - (*m_Nodes[0])[0];
  const PetscReal dy = (*m_Nodes[1])[1] - (*m_Nodes[0])[1];

  const PetscReal nx = dx / L;
  const PetscReal ny = dy / L;

  // --- create transformationmatrix
  T(0, 0) = nx;
  T(0, 1) = ny;
  T(1, 0) = -ny;
  T(1, 1) = nx;
  T(2, 2) = 1.;
  T(3, 3) = nx;
  T(3, 4) = ny;
  T(4, 3) = -ny;
  T(4, 4) = nx;
  T(5, 5) = 1.;

  if (transpose == true) {
    std::swap(T(0, 1), T(1, 0));
    std::swap(T(3, 4), T(4, 3));
  }
}

std::vector<eKnownDofs> cElementStructureBeam::getDofs(void) const {
  std::vector<eKnownDofs> res(3);

  res[0] = disp_x1;
  res[1] = disp_x2;
  res[2] = disp_w3;

  return res;
}

PetscReal cElementStructureBeam::getPsi(void) const {
  if (m_UseTheory == Bernoulli) {
    return 1.0;
  } else {
    const PetscScalar E = m_Material->getEOmega();
    const PetscReal k = m_Material->getKs();
    const PetscReal I = m_Material->getI();
    const PetscReal L = getLength();
    const PetscReal nu = m_Material->getNu();
    const PetscScalar G = E / (2. * (1. + nu));
    const PetscReal A = m_Material->getA();

#ifdef PETSC_USE_COMPLEX
    const PetscReal res =
        1. / (1. + 12. * E.real() * I / (L * L * G.real() * k * A));
#else
    const PetscReal res = 1. / (1. + 12. * E * I / (L * L * G * k * A));
#endif

    return res;
  }
}

void cElementStructureBeam::assembleStiffnessMatrix(cElementMatrix &KM,
                                                    Vec *x = NULL,
                                                    Vec *dx = NULL) {
  cElementMatrix K(6, 6);

  const PetscScalar E = m_Material->getEOmega();
  const PetscReal A = m_Material->getA();
  const PetscReal I = m_Material->getI();
  const PetscReal L = getLength();
  const PetscScalar EA = E * A;
  const PetscScalar EI = E * I;
  const PetscReal PSI = getPsi();

  // -------------------------------------------------------------------------
  //   stiffness matrix: combination of a rod and a Timoshenko beam element
  K(0, 0) = EA / L;
  K(0, 3) = -EA / L;

  K(1, 1) = (EI / (L * L * L)) * 12. * PSI;
  K(1, 2) = (EI / (L * L * L)) * 6. * L * PSI;
  K(1, 4) = (EI / (L * L * L)) * -12. * PSI;
  K(1, 5) = (EI / (L * L * L)) * 6. * L * PSI;

  K(2, 1) = (EI / (L * L * L)) * 6. * L * PSI;
  K(2, 2) = (EI / (L * L * L)) * L * L * (1. + 3. * PSI);
  K(2, 4) = (EI / (L * L * L)) * -6. * L * PSI;
  K(2, 5) = (EI / (L * L * L)) * L * L * (-1. + 3. * PSI);

  K(3, 0) = -EA / L;
  K(3, 3) = EA / L;

  K(4, 1) = (EI / (L * L * L)) * -12. * PSI;
  K(4, 2) = (EI / (L * L * L)) * -6. * L * PSI;
  K(4, 4) = (EI / (L * L * L)) * 12. * PSI;
  K(4, 5) = (EI / (L * L * L)) * -6. * L * PSI;

  K(5, 1) = (EI / (L * L * L)) * 6. * L * PSI;
  K(5, 2) = (EI / (L * L * L)) * L * L * (-1. + 3. * PSI);
  K(5, 4) = (EI / (L * L * L)) * -6. * L * PSI;
  K(5, 5) = (EI / (L * L * L)) * L * L * (1. + 3. * PSI);

  // --- compute KM = T^t * K * T
  cMatrix T(6, 6);
  computeTransformationMatrix(T);
  BT_C_B_wdJ(KM, T, K, 1.);
}

void cElementStructureBeam::assembleMassMatrix(cElementMatrix &MM) {
  cElementMatrix M(6, 6);

  const PetscReal A = m_Material->getA();
  const PetscReal I = m_Material->getI();
  const PetscReal rho = m_Material->getRho();
  const PetscReal L = getLength();
  const PetscReal PSI = getPsi();

  // --- mass matrix (rod and Timoshenko beam element)
  M(0, 0) = rho * A * L / 3.;
  M(0, 3) = rho * A * L / 6.;

  M(1, 1) = (rho * A * L / 840.) * (280. + 28. * PSI + 4. * PSI * PSI);
  M(1, 2) = (rho * A * L / 840.) * (L * (35. + 7. * PSI + 2. * PSI * PSI));
  M(1, 4) = (rho * A * L / 840.) * (140. - 28. * PSI - 4. * PSI * PSI);
  M(1, 5) = (rho * A * L / 840.) * (-L * (35. - 7. * PSI - 2. * PSI * PSI));

  M(2, 1) = (rho * A * L / 840.) * (L * (35. + 7. * PSI + 2. * PSI * PSI));
  M(2, 2) = (rho * A * L / 840.) * (L * L * (7. + PSI * PSI));
  M(2, 4) = (rho * A * L / 840.) * (L * (35. - 7. * PSI - 2. * PSI * PSI));
  M(2, 5) = (rho * A * L / 840.) * (-L * L * (7. - PSI * PSI));

  M(3, 0) = rho * A * L / 6.;
  M(3, 3) = rho * A * L / 3.;

  M(4, 1) = (rho * A * L / 840.) * (140. - 28. * PSI - 4. * PSI * PSI);
  M(4, 2) = (rho * A * L / 840.) * (L * (35. - 7. * PSI - 2. * PSI * PSI));
  M(4, 4) = (rho * A * L / 840.) * (280. + 28. * PSI + 4. * PSI * PSI);
  M(4, 5) = (rho * A * L / 840.) * (-L * (35. + 7. * PSI + 2. * PSI * PSI));

  M(5, 1) = (rho * A * L / 840.) * (-L * (35. - 7. * PSI - 2. * PSI * PSI));
  M(5, 2) = (rho * A * L / 840.) * (-L * L * (7. - PSI * PSI));
  M(5, 4) = (rho * A * L / 840.) * (-L * (35. + 7. * PSI + 2. * PSI * PSI));
  M(5, 5) = (rho * A * L / 840.) * (L * L * (7. + PSI * PSI));

  M(1, 1) += (rho * I / (30. * L)) * (36. * PSI * PSI);
  M(1, 2) += (rho * I / (30. * L)) * (-L * (15. * PSI - 18. * PSI * PSI));
  M(1, 4) += (rho * I / (30. * L)) * (-36. * PSI * PSI);
  M(1, 5) += (rho * I / (30. * L)) * (-L * (15. * PSI - 18. * PSI * PSI));

  M(2, 1) += (rho * I / (30. * L)) * (-L * (15. * PSI - 18. * PSI * PSI));
  M(2, 2) +=
      (rho * I / (30. * L)) * (L * L * (10. - 15. * PSI + 9. * PSI * PSI));
  M(2, 4) += (rho * I / (30. * L)) * (L * (15. * PSI - 18. * PSI * PSI));
  M(2, 5) +=
      (rho * I / (30. * L)) * (L * L * (5. - 15. * PSI + 9. * PSI * PSI));

  M(4, 1) += (rho * I / (30. * L)) * (-36. * PSI * PSI);
  M(4, 2) += (rho * I / (30. * L)) * (L * (15. * PSI - 18. * PSI * PSI));
  M(4, 4) += (rho * I / (30. * L)) * (36. * PSI * PSI);
  M(4, 5) += (rho * I / (30. * L)) * (L * (15. * PSI - 18. * PSI * PSI));

  M(5, 1) += (rho * I / (30. * L)) * (-L * (15. * PSI - 18. * PSI * PSI));
  M(5, 2) +=
      (rho * I / (30. * L)) * (L * L * (5. - 15. * PSI + 9. * PSI * PSI));
  M(5, 4) += (rho * I / (30. * L)) * (L * (15. * PSI - 18. * PSI * PSI));
  M(5, 5) +=
      (rho * I / (30. * L)) * (L * L * (10. - 15. * PSI + 9. * PSI * PSI));

  // --- 'exact' mass matrix of an Euler Bernoulli beam
  //   const PetscScalar factor = rho * A * L / 420.;
  //   M(0,0) = 140.;
  //   M(0,3) =  70.;
  //   M(3,0) =  70.;
  //   M(3,3) = 140.;
  //   M(1,1) = 156.;
  //   M(1,2) =  22. * L;
  //   M(1,4) =  54.;
  //   M(1,5) = -13. * L;
  //   M(2,1) =  22. * L;
  //   M(2,2) =   4. * L * L;
  //   M(2,4) =  13. * L;
  //   M(2,5) =  -3. * L * L;
  //   M(4,1) =  54.;
  //   M(4,2) =  13. * L;
  //   M(4,4) = 156.;
  //   M(4,5) = -22. * L;
  //   M(5,1) = -13. * L;
  //   M(5,2) =  -3. * L * L;
  //   M(5,4) = -22. * L;
  //   M(5,5) =   4. * L * L;

  //   infam::scale(M, factor);

  // --- compute MM = T^t * M * T
  cMatrix T(6, 6);
  computeTransformationMatrix(T);
  BT_C_B_wdJ(MM, T, M, 1.);
}

void cElementStructureBeam::assembleLoadVector(cElementVector &LV,
                                               cElementMatrix &KM,
                                               Vec *x = NULL, Vec *dx = NULL) {
  // --- leave if no load is applied
  //     to this element
  if (m_ElementLoads.size() != 0) {
    // --- Px and Py are constant elementloads whose x and y direction
    //     correspond to the GLOBAL coordinate system
    PetscScalar Px = 0.;
    PetscScalar Py = 0.;
    const PetscReal L = getLength();

    // --- first the transformation matrix T is used
    //     to compute the local forces of this element
    cMatrix T(6, 6);
    computeTransformationMatrix(T);

    // --- load vector in LOCAL coordinate system
    cElementVector LVL(6);

    std::multimap<short, cElementLoad *>::iterator it;
    cElementLoadStructure *ptrELS = NULL;

    for (it = m_ElementLoads.begin(); it != m_ElementLoads.end(); it++) {
      ptrELS = dynamic_cast<cElementLoadStructure *>(it->second);

      if (ptrELS != NULL) {
        // --- we assume the load to be constant
        Px = ptrELS->getForceComponent(m_Nodes[0]->getId(), 0);
        Py = ptrELS->getForceComponent(m_Nodes[0]->getId(), 1);

        // --- we need the force components in local
        //     coordinate directions
        const PetscScalar PLx = T(0, 0) * Px + T(0, 1) * Py;
        const PetscScalar PLy = T(1, 0) * Px + T(1, 1) * Py;

        LVL[0] += 1. / 2. * PLx * L;
        LVL[1] += 1. / 2. * PLy * L;
        LVL[2] += -1. / 12. * PLy * L * L;
        LVL[3] += 1. / 2. * PLx * L;
        LVL[4] += 1. / 2. * PLy * L;
        LVL[5] += 1. / 12. * PLy * L * L;
      }
      ptrELS = NULL;
    }

    // --- setup transformationmatrix
    //     and compute LV += T^t * LVL
    computeTransformationMatrix(T, true);  // now we need T^t
    infam::mult(T, LVL, LV);
  }
}

std::ostream &cElementStructureBeam::write(std::ostream &os) const {
  os << "beam element (Euler-Bernoulli)" << std::endl;
  os << "  no. nodes : " << getNumberOfNodes();
  os << "  dofs/node : " << getNumberOfDofsPerNode();
  for (int k = 0; k < getNumberOfNodes(); k++) os << *(m_Nodes[k]);

  return os;
}

std::ostream &cElementStructureBeam::writeXml(std::ostream &os) const {
  os << "<Beam>";
  os << "<Id>" << getId() << "</Id>";
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</Beam>";
  return os;
}
