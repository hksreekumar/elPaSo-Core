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

#include "elementstructurebeam3d.h"

#include <iomanip>


cElementStructureBeam3D::cElementStructureBeam3D(short NumberOfNodes,
                                                 short NumberOfDofsPerNode,
                                                 short NumberOfGaussPoints,
                                                 const eUseBeamTheory &Theory)
    : cElementStructureLinear(NumberOfNodes, NumberOfDofsPerNode,
                              NumberOfGaussPoints) {
  m_UseTheory = Theory;
}

cElementStructureBeam3D::cElementStructureBeam3D(
    const cElementStructureBeam3D &other)
    : cElementStructureLinear(other) {
  m_UseTheory = other.getBeamTheory();
}

cElementStructureBeam3D::~cElementStructureBeam3D() {
  // empty
}

cVector cElementStructureBeam3D::getGlobalNormalVector(void) const {
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

  const PetscReal L = getLength();
  PetscReal nx = (new_b[0] - new_a[0]) / L;
  PetscReal ny = (new_b[1] - new_a[1]) / L;
  PetscReal nz = (new_b[2] - new_a[2]) / L;

  //// --- local normal vector n^T = [ 0, 1, 0 ];
  ////     globalNormal = T^t * n
  cVector globalNormal(3);

  globalNormal[0] = ny;  // local nx
  globalNormal[1] = nx;  // local ny
  globalNormal[2] = nz;  // local nz

  // infam::scale( globalNormal, 1. / globalNormal.abs() );
  return globalNormal;
}

PetscReal cElementStructureBeam3D::getPsi(void) const {
  if (m_UseTheory == Bernoulli) {
    return 1.0;
  } else {
    const PetscScalar E = m_Material->getEOmega();
    const PetscReal k = m_Material->getKs();
    const PetscReal L = getLength();
    const PetscReal nu = m_Material->getNu();
    const PetscScalar G = E / (2. * (1. + nu));
    const PetscReal A = m_Material->getA();

#ifdef PETSC_USE_COMPLEX
    const PetscReal res = 12. * E.real() / (L * L * G.real() * k * A);
#else
    const PetscReal res = 12. * E / (L * L * G * k * A);
#endif
    return res;
  }
}
