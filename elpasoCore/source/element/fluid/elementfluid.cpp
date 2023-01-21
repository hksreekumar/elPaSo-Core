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

#include "elementfluid.h"

cElementFluid::cElementFluid(short NumberOfNodes, short NumberOfGaussPoints)
    : cElementFEM(NumberOfNodes, 1, NumberOfGaussPoints) {
  m_Material = NULL;
  m_nl = true;
}

cElementFluid::cElementFluid(const cElementFluid &other) : cElementFEM(other) {
  setMaterial(other.getMaterial());
  m_nl = other.isNonlinearElement();
}

cElementFluid::~cElementFluid() { m_Material = NULL; }

cMaterial *cElementFluid::getMaterial(void) const { return m_Material; }

std::vector<eKnownDofs> cElementFluid::getDofs(void) const {
  std::vector<eKnownDofs> res(1);
  res[0] = fluid;

  return res;
}

void cElementFluid::setMaterial(cMaterial *ptrMaterial) {
  if (ptrMaterial == NULL)
    throw cException("ptrMaterial is NULL", __FILE__, __LINE__);

  m_Material = dynamic_cast<cMaterialFluid *>(ptrMaterial);

  if (m_Material == NULL) {
    message(
        "  wrong material type (FLUID element no. %d, id of material : %d\n",
        getId(), ptrMaterial->getId());
    throw cException("invalid materialset", __FILE__, __LINE__);
  }
}

void cElementFluid::assembleDynamicStiffnessMatrix(const PetscReal &omega,
                                                   cElementMatrix &EM) {
  const int nnod = getNumberOfNodes();
  cElementMatrix KK(nnod, nnod);
  cElementMatrix MM(nnod, nnod);

  assembleMassMatrix(MM);
  infam::scale(MM, -omega * omega);
  // std::cout<<"ACHTUNG cElementFluid::assembleDynamicStiffnessMatrix vector
  // a,b nicht definiert!!\n";
  assembleStiffnessMatrix(KK, NULL, NULL);
  evaluateImpedance(omega, KK);

  infam::add(MM, EM);
  infam::add(KK, EM);

  // --- divide by rho*omega^2 for symmetric formulation
  const PetscScalar factor = 1. / (m_Material->getRhoOmega() * omega * omega);
  infam::scale(EM, factor);
}
