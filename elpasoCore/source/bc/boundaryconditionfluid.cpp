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

#include "boundaryconditionfluid.h"

/*BEGIN_NO_COVERAGE*/
cBoundaryConditionFluid::cBoundaryConditionFluid() { setPressure(0.); }

cBoundaryConditionFluid::cBoundaryConditionFluid(
    const cBoundaryConditionFluid& other)
    : cBoundaryCondition(other) {
  setPressure(other.getPressure());
}

cBoundaryConditionFluid::~cBoundaryConditionFluid() {
  // empty
}
/*END_NO_COVERAGE*/

bool cBoundaryConditionFluid::checkIfFixed(eKnownDofs dof) const {
  if (dof == fluid)
    return true;
  else
    return false;
}

bool cBoundaryConditionFluid::checkIfFixed(int dof) const {
  if (dof == 16)
    return true;
  else
    return false;
}

PetscScalar cBoundaryConditionFluid::getPrescribedValue(eKnownDofs dof,
                                                        cPoint* curPt,
                                                        PetscReal omega) const {
  if (dof == fluid)
    return (PetscScalar)getPressure();
  else
    return 0.;
}

/*BEGIN_NO_COVERAGE*/
std::ostream& cBoundaryConditionFluid::write(std::ostream& os) const {
  os << "Nr. ";
  os.width(2);
  os << getId();
  os << "  (Fluid-Randbedingung)" << std::endl;
  os << std::setw(5) << getPressure() << std::endl;
  return os;
}

std::istream& cBoundaryConditionFluid::read(std::istream& is) {
  cId::read(is);
  is >> m_Identifier;
  is >> m_Pressure;

  return is;
}

std::ostream& cBoundaryConditionFluid::writeXml(std::ostream& os) const {
  os << "<NodeBC Type=\"fluid\" Name=\"" << getIdentifier() << "\">"
     << std::endl;
  os << "  <Id>" << getId() << "</Id>" << std::endl;
  os << "  <P>" << getPressure() << "</P>" << std::endl;
  os << "</NodeBC>" << std::endl;

  return os;
}
/*END_NO_COVERAGE*/