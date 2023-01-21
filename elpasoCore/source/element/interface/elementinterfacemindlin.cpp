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

#include "elementinterfacemindlin.h"

cElementInterfaceMindlin::cElementInterfaceMindlin(short NumberOfNodes)
    : cElementInterfaceStructural(NumberOfNodes) {
  m_MaterialFluid = NULL;
}

cElementInterfaceMindlin::cElementInterfaceMindlin(
    const cElementInterfaceMindlin &other)
    : cElementInterfaceStructural(other) {
  setMaterial(other.getMaterial(0), 0);
}

cElementInterfaceMindlin::~cElementInterfaceMindlin() {
  m_MaterialFluid = NULL;
}

void cElementInterfaceMindlin::setMaterial(cMaterial *ptrMaterial,
                                           const int &domain) {
  if (domain == 0) {
    m_MaterialFluid = dynamic_cast<cMaterialFluid *>(ptrMaterial);
    if (m_MaterialFluid == NULL) {
      PetscPrintf(PETSC_COMM_SELF,
                  "ERROR: invalid material %d for interfaceelement %d\n",
                  ptrMaterial->getId(), getId());
      ExitApp();
    }
  }
}

cMaterial *cElementInterfaceMindlin::getMaterial(const int &domain) const {
  if (domain == 0)
    return m_MaterialFluid;
  else
    return NULL;
}

eElementShape cElementInterfaceMindlin::getElementShape(void) const {
  if ((getNumberOfNodes() == 4) || (getNumberOfNodes() == 9))
    return Quadrilateral;
  else {
    PetscPrintf(PETSC_COMM_SELF,
                "ERROR: No interface element for %d nodes implemented",
                getNumberOfNodes());
    ExitApp();
  }

  return Quadrilateral;
}

std::vector<eKnownDofs> cElementInterfaceMindlin::getDofsStructure(void) const {
  std::vector<eKnownDofs> res(3);

  res[0] = disp_x1;
  res[1] = disp_x2;
  res[2] = disp_x3;

  return res;
}

void cElementInterfaceMindlin::computeCouplingMatrix(cElementMatrix &C) {
  const PetscInt nnod = getNumberOfNodes();  // number of nodes of this element
  const int ngp = (nnod == 4 ? 2 : 3);       // number of Gauss points used
  cPoint gp;                                 // current Gauss point
  PetscReal N[9], Nxi[9], Neta[9];           // evaluated testfunctions
  cVector a1(3);                             // basis of tangential plane at ...
  cVector a2(3);                             // ... current Gauss point
  PetscReal detJac;                          // Jacobian
  PetscReal weight;                          // weights of numerical integration

  // -------------------------------------------------------------------------
  //   perform numerical integration
  // -------------------------------------------------------------------------
  for (int n = 0; n < ngp * ngp; n++) {
    // --- evaluate testfunctions
    gp = m_GaussPoints.getGaussPoint2D(ngp, n);
    weight = m_GaussPoints.getGaussWeight2D(ngp, n);

    for (int k = 0; k < nnod; k++) {
      N[k] = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k,
                                                    N_fun, gp);
      Nxi[k] = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k,
                                                      N_xi, gp);
      Neta[k] = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k,
                                                       N_eta, gp);
    }

    // --- compute Jacobian
    a1.setValue(0.0);
    a2.setValue(0.0);

    for (int k = 0; k < nnod; k++) {
      for (int d = 0; d < 3; d++) {
        a1[d] += Nxi[k] * (*m_Nodes[k])[d];
        a2[d] += Neta[k] * (*m_Nodes[k])[d];
      }
    }

    detJac = sqrt(a1.abs2() * a2.abs2() - (a1.dot(a2)) * (a1.dot(a2)));

    for (int z = 0; z < nnod; z++)
      for (int s = 0; s < nnod; s++)
        C(3 * z + 2, s) += N[z] * N[s] * weight * detJac;
  }

  // --- transformation into global coordinate system
  //     this function call is only needed when the structure is discretized
  //     by shell elements.
  if (m_MatchingNodes[0]->checkIfActive(disp_x1) == true) transformC(C);
}

void cElementInterfaceMindlin::transformC(cElementMatrix &C) {
  const PetscInt nnod = getNumberOfNodes();  // number of nodes of this element
  cPoint gp;                                 // current Gauss point

  // --- coordinates of element's nodes
  PetscReal m_Xi[9] = {-1., +1., +1., -1., 0., +1., 0., -1., 0.};
  PetscReal m_Eta[9] = {-1., -1., +1., +1., -1., 0., +1., 0., 0.};
  cVector normal(3);          // normal vector at one element node
  PetscReal Nxi[9], Neta[9];  // evaluated testfunctions
  cVector a1(3);              // basis of tangential plane at ...
  cVector a2(3);              // ... current Gauss point

  for (int n = 0; n < nnod; n++) {
    gp[0] = m_Xi[n];
    gp[1] = m_Eta[n];
    gp[2] = 0.;

    for (int k = 0; k < nnod; k++) {
      Nxi[k] = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k,
                                                      N_xi, gp);
      Neta[k] = m_ShapeFunctions.evaluateShapeFunction(Quadrilateral, nnod, k,
                                                       N_eta, gp);
    }

    // --- compute normal vector at one elementnode:
    //     crossproduct of tangential vectors
    a1.setValue(0.0);
    a2.setValue(0.0);

    for (int k = 0; k < nnod; k++) {
      for (int d = 0; d < 3; d++) {
        a1[d] += Nxi[k] * (*m_Nodes[k])[d];
        a2[d] += Neta[k] * (*m_Nodes[k])[d];
      }
    }

    infam::cross_product(a1, a2,
                         normal);  // compute crossproduct: normal = a1 x a2
    infam::scale(normal, 1. / normal.abs());  // scale vector to abs()=1

    // compute the projection of the normal derivative on the global
    // coordinate directions
    //  e.g. e_3 \cdot n = [0 0 1] \cdot [n_1 n_2 n_3] = n_3
    for (int s = 0; s < nnod; s++) {
      C(3 * n, s) = normal[0] * C(3 * n + 2, s);
      C(3 * n + 1, s) = normal[1] * C(3 * n + 2, s);
      C(3 * n + 2, s) =
          normal[2] * C(3 * n + 2, s);  // now the value can be modified
    }
  }
}

std::ostream &cElementInterfaceMindlin::write(std::ostream &os) const {
  os << "Interface Helmholtz <-> Mindlin" << std::endl;
  cElementInterface::write(os);
  return os;
}

std::ostream &cElementInterfaceMindlin::writeXml(std::ostream &os) const {
  os << "<InterfaceMindlin n=\"" << getNumberOfNodes() << "\">";
  os << std::endl;
  os << "<Id>" << getId() << "</Id>";
  os << "<MatF>" << getMaterial(0)->getId() << "</MatF>";
  os << "<ori>" << std::showpoint << getOrientation() << "</ori>" << std::endl;
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<NF>" << m_Nodes[k]->getId() << "</NF>";
  os << std::endl;
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<NS>" << getMatchingNode(k)->getId() << "</NS>";
  os << std::endl;
  os << "</InterfaceMindlin>";
  return os;
}
