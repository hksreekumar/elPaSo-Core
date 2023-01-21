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

#include "elementinterfacekirch.h"

cElementInterfaceKirchhoff::cElementInterfaceKirchhoff(short NumberOfNodes)
    : cElementInterfaceStructural(NumberOfNodes) {
  m_MaterialFluid = NULL;
}

cElementInterfaceKirchhoff::cElementInterfaceKirchhoff(
    const cElementInterfaceKirchhoff &other)
    : cElementInterfaceStructural(other) {
  setMaterial(other.getMaterial(0), 0);
}

cElementInterfaceKirchhoff::~cElementInterfaceKirchhoff() {
  m_MaterialFluid = NULL;
}

void cElementInterfaceKirchhoff::setMaterial(cMaterial *ptrMaterial,
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

cMaterial *cElementInterfaceKirchhoff::getMaterial(const int &domain) const {
  if (domain == 0)
    return m_MaterialFluid;
  else
    return NULL;
}

eElementShape cElementInterfaceKirchhoff::getElementShape(void) const {
  if (getNumberOfNodes() == 4) return Quadrilateral;

  return Quadrilateral;
}

std::vector<eKnownDofs> cElementInterfaceKirchhoff::getDofsStructure(
    void) const {
  std::vector<eKnownDofs> res(4);

  res[0] = disp_x3;
  res[1] = disp_dwdx;
  res[2] = disp_dwdy;
  res[3] = disp_dwdxy;

  return res;
}

void cElementInterfaceKirchhoff::findLowerCorner(void) {
  std::vector<double>::iterator itMinX, itMinY;
  std::vector<double> x(4);
  std::vector<double> y(4);

  for (int k = 0; k < 4; k++) {
    x[k] = (*m_Nodes[k])[0];
    y[k] = (*m_Nodes[k])[1];
  }

  itMinX = std::min_element(x.begin(), x.end());
  itMinY = std::min_element(y.begin(), y.end());

  // --- find the node that contains both minimum values
  int newFirst = 0;
  for (int k = 0; k < getNumberOfNodes(); k++) {
    if ((std::abs((*m_Nodes[k])[0] - (*itMinX)) < cstGeomEps) &&
        (std::abs((*m_Nodes[k])[1] - (*itMinY)) < cstGeomEps)) {
      newFirst = k;
    }
  }

  // --- rotate the element
  if (newFirst != 0) {
    std::rotate(m_Nodes.begin(), m_Nodes.begin() + newFirst, m_Nodes.end());
    std::rotate(m_MatchingNodes.begin(), m_MatchingNodes.begin() + newFirst,
                m_MatchingNodes.end());
  }

  // --- okay. Now the first node of the interface element
  //     is the lower left one. But if the normals of fluid elements face
  //     and the structural element point into opposite direction, we
  //     can not compute lx and ly. Therefore we change the order of the
  //     nodes of the element.
  if (getOrientation() > 0.) {
    std::swap(m_Nodes[1], m_Nodes[3]);
    std::swap(m_MatchingNodes[1], m_MatchingNodes[3]);
  }
}

void cElementInterfaceKirchhoff::computeCouplingMatrix(cElementMatrix &C) {
  PetscReal a =
      (*m_Nodes[1])[0] - (*m_Nodes[0])[0];  // length of element x-direction
  PetscReal b =
      (*m_Nodes[3])[1] - (*m_Nodes[0])[1];  // length of element y-direction

  if ((a < cstGeomEps) || (b < cstGeomEps)) {
    a = (*m_Nodes[3])[0] - (*m_Nodes[0])[0];  // length of element x-direction
    b = (*m_Nodes[1])[1] - (*m_Nodes[0])[1];  // length of element y-direction

    //   message("  Error in interface element no. %d\n", getId());
    /*    message("    a = %f\n", a);
       message("    b = %f\n", b);
      ExitApp();*/
  }

  C(0, 0) = 0.49e2 / 0.400e3 * a * b;
  C(0, 1) = 0.21e2 / 0.400e3 * a * b;
  C(0, 2) = 0.9e1 / 0.400e3 * a * b;
  C(0, 3) = 0.21e2 / 0.400e3 * a * b;
  C(1, 0) = 0.7e1 / 0.400e3 * a * a * b;
  C(1, 1) = 0.7e1 / 0.600e3 * a * a * b;
  C(1, 2) = a * a * b / 0.200e3;
  C(1, 3) = 0.3e1 / 0.400e3 * a * a * b;
  C(2, 0) = 0.7e1 / 0.400e3 * a * b * b;
  C(2, 1) = 0.3e1 / 0.400e3 * a * b * b;
  C(2, 2) = a * b * b / 0.200e3;
  C(2, 3) = 0.7e1 / 0.600e3 * a * b * b;
  C(3, 0) = a * a * b * b / 0.400e3;
  C(3, 1) = a * a * b * b / 0.600e3;
  C(3, 2) = a * a * b * b / 0.900e3;
  C(3, 3) = a * a * b * b / 0.600e3;
  C(4, 0) = 0.21e2 / 0.400e3 * a * b;
  C(4, 1) = 0.49e2 / 0.400e3 * a * b;
  C(4, 2) = 0.21e2 / 0.400e3 * a * b;
  C(4, 3) = 0.9e1 / 0.400e3 * a * b;
  C(5, 0) = -0.7e1 / 0.600e3 * a * a * b;
  C(5, 1) = -0.7e1 / 0.400e3 * a * a * b;
  C(5, 2) = -0.3e1 / 0.400e3 * a * a * b;
  C(5, 3) = -a * a * b / 0.200e3;
  C(6, 0) = 0.3e1 / 0.400e3 * a * b * b;
  C(6, 1) = 0.7e1 / 0.400e3 * a * b * b;
  C(6, 2) = 0.7e1 / 0.600e3 * a * b * b;
  C(6, 3) = a * b * b / 0.200e3;
  C(7, 0) = -a * a * b * b / 0.600e3;
  C(7, 1) = -a * a * b * b / 0.400e3;
  C(7, 2) = -a * a * b * b / 0.600e3;
  C(7, 3) = -a * a * b * b / 0.900e3;
  C(8, 0) = 0.9e1 / 0.400e3 * a * b;
  C(8, 1) = 0.21e2 / 0.400e3 * a * b;
  C(8, 2) = 0.49e2 / 0.400e3 * a * b;
  C(8, 3) = 0.21e2 / 0.400e3 * a * b;
  C(9, 0) = -a * a * b / 0.200e3;
  C(9, 1) = -0.3e1 / 0.400e3 * a * a * b;
  C(9, 2) = -0.7e1 / 0.400e3 * a * a * b;
  C(9, 3) = -0.7e1 / 0.600e3 * a * a * b;
  C(10, 0) = -a * b * b / 0.200e3;
  C(10, 1) = -0.7e1 / 0.600e3 * a * b * b;
  C(10, 2) = -0.7e1 / 0.400e3 * a * b * b;
  C(10, 3) = -0.3e1 / 0.400e3 * a * b * b;
  C(11, 0) = a * a * b * b / 0.900e3;
  C(11, 1) = a * a * b * b / 0.600e3;
  C(11, 2) = a * a * b * b / 0.400e3;
  C(11, 3) = a * a * b * b / 0.600e3;
  C(12, 0) = 0.21e2 / 0.400e3 * a * b;
  C(12, 1) = 0.9e1 / 0.400e3 * a * b;
  C(12, 2) = 0.21e2 / 0.400e3 * a * b;
  C(12, 3) = 0.49e2 / 0.400e3 * a * b;
  C(13, 0) = 0.3e1 / 0.400e3 * a * a * b;
  C(13, 1) = a * a * b / 0.200e3;
  C(13, 2) = 0.7e1 / 0.600e3 * a * a * b;
  C(13, 3) = 0.7e1 / 0.400e3 * a * a * b;
  C(14, 0) = -0.7e1 / 0.600e3 * a * b * b;
  C(14, 1) = -a * b * b / 0.200e3;
  C(14, 2) = -0.3e1 / 0.400e3 * a * b * b;
  C(14, 3) = -0.7e1 / 0.400e3 * a * b * b;
  C(15, 0) = -a * a * b * b / 0.600e3;
  C(15, 1) = -a * a * b * b / 0.900e3;
  C(15, 2) = -a * a * b * b / 0.600e3;
  C(15, 3) = -a * a * b * b / 0.400e3;
}

std::ostream &cElementInterfaceKirchhoff::write(std::ostream &os) const {
  os << "Interface Helmholtz <-> Kirchhoff" << std::endl;
  cElementInterface::write(os);
  return os;
}

std::ostream &cElementInterfaceKirchhoff::writeXml(std::ostream &os) const {
  os << "<InterfaceKirch n=\"" << getNumberOfNodes() << "\">";
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
  os << "</InterfaceKirch>";
  return os;
}
