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

#include "boundaryconditionfactory.h"

#include "../misc/parser/femparserinterface.h"
#include "boundaryconditionfluid.h"
#include "boundaryconditionstructure.h"

cBoundaryConditionFactory::cBoundaryConditionFactory() {
  // empty
}

cBoundaryConditionFactory::~cBoundaryConditionFactory() {
  // empty
}

cBoundaryCondition* cBoundaryConditionFactory::createBoundaryCondition(
    std::string _bcIdentifier, int _id, cFemParserInterface* _parser) {
  cBoundaryCondition* ptrBoundCond = 0;
  std::stringstream data;

  if (_bcIdentifier == "structure") {
    ptrBoundCond = new cBoundaryConditionStructure;
    ParserNodeConstraintStructureData readData =
        _parser->getNodeConstraintsStructureData(_id);

    data << readData.nconstraintId << " " << readData.nconstraintName << " ";
    for (int k = 0; k < readData.tagcount; k++)
      data << readData.nconstraintFlag[k] << " " << readData.nconstraintVals[k]
           << " ";
    data >> *ptrBoundCond;
  } else if (_bcIdentifier == "acoustic") {
    ptrBoundCond = new cBoundaryConditionFluid;
    ParserNodeConstraintAcousticData readData =
        _parser->getNodeConstraintsAcousticData(_id);

    data << readData.nconstraintId << " " << readData.nconstraintName << " "
         << readData.nconstraintValue << " ";
    data >> *ptrBoundCond;
  }

  return ptrBoundCond;
}
