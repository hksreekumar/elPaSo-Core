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

#include "elementinterfacefactory.h"

#include "../../misc/parser/femparserinterface.h"
#include "elementinterfacekirch.h"
#include "elementinterfacemindlin.h"

cElementInterfaceFactory::cElementInterfaceFactory() {
  // empty
}

cElementInterfaceFactory::~cElementInterfaceFactory() {
  // empty
}

bool cElementInterfaceFactory::createElementInterface(
    std::string Type, AppInterfaceElements _data, void* userData) {
  // Identify element type and direct to element definition
  if (Type == "InterfaceKirch")
    sElementInterfaceKirch(&_data, userData);
  else if (Type == "InterfaceMindlin")
    sElementInterfaceMindlin(&_data, userData);
  else
    return false;

  return true;
}

void cElementInterfaceFactory::sElementInterfaceKirch(AppInterfaceElements* App,
                                                      void* userData) {
  cElementInterfaceKirchhoff* ptr = new cElementInterfaceKirchhoff(App->nnod);
  parseSingleInterfaceElement(App, ptr, userData);
  ptr->findLowerCorner();
  ptr = NULL;
}

void cElementInterfaceFactory::sElementInterfaceMindlin(
    AppInterfaceElements* App, void* userData) {
  cElementInterface* ptr = new cElementInterfaceMindlin(App->nnod);
  parseSingleInterfaceElement(App, ptr, userData);
  ptr = NULL;
}

void cElementInterfaceFactory::parseSingleInterfaceElement(
    AppInterfaceElements* App, cElementInterface* ptrElement, void* userData) {
  // -------------------------------------------------------------------------
  //   WARNING: the sequence of the entries is important!
  //            first the Id, followed by MaterialId and the nodes
  // -------------------------------------------------------------------------
  int ElemId, NodeId, MatF, StrId, FluId;
  double ori;
  try {
    MatF = App->MatF;
    ElemId = App->dataspace[0];
    ori = App->dataspace[1];

    // -------------------------------------------------------------------------
    //   construct element of data read
    // -------------------------------------------------------------------------
    cMesh* MyMesh = ((cProblem*)userData)->getMesh();
    ptrElement->setId(ElemId);
    ptrElement->setOrientation(ori);

    ptrElement->setMaterial(MyMesh->getMaterial(MatF), 0);

    for (int i = 0; i < ptrElement->getNumberOfNodes(); i++) {
      NodeId = App->dataspace[2 + ptrElement->getNumberOfNodes() +
                              i];  // Fluid nodes
      ptrElement->setNode(i, MyMesh->getNode(NodeId));

      NodeId = App->dataspace[2 + i];  // Structural nodes
      ptrElement->setMatchingNode(i, MyMesh->getNode(NodeId));
    }
    MyMesh->insertInterfaceElement(ptrElement);
  } catch (...) {
    PetscPrintf(PETSC_COMM_SELF, "unknown error in file %s, line %s catched.",
                __FILE__, __LINE__);
    ExitApp();
  }
}
