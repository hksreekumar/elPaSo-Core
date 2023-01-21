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

#include "boundaryconditionstructure.h"

/*BEGIN_NO_COVERAGE*/
cBoundaryConditionStructure::cBoundaryConditionStructure()
    : m_Value(eNumberOfDofs, 0.0) {
  m_Fixed.reset();
}

cBoundaryConditionStructure::cBoundaryConditionStructure(
    const cBoundaryConditionStructure& other)
    : cBoundaryCondition(other),
      m_Fixed(other.getFixedArray()),
      m_Value(other.getValueArray()) {
  // empty
}

cBoundaryConditionStructure::~cBoundaryConditionStructure() {
  // empty
}
/*END_NO_COVERAGE*/

bool cBoundaryConditionStructure::checkIfFixed(eKnownDofs dof) const {
  if (dof == disp_x1) return m_Fixed.test(0);
  if (dof == disp_x2) return m_Fixed.test(1);
  if (dof == disp_x3) return m_Fixed.test(2);
  if (dof == disp_w1) return m_Fixed.test(3);
  if (dof == disp_w2) return m_Fixed.test(4);
  if (dof == disp_w3) return m_Fixed.test(5);
  if (dof == pore0) return m_Fixed.test(6);
  if (dof == pore1) return m_Fixed.test(7);
  if (dof == pore2) return m_Fixed.test(8);
  if (dof == pore3) return m_Fixed.test(9);
  if (dof == disp_xd3) return m_Fixed.test(10);
  if (dof == disp_wd1) return m_Fixed.test(11);
  if (dof == disp_wd2) return m_Fixed.test(12);
  if (dof == disp_dwdx) return m_Fixed.test(13);
  if (dof == disp_dwdy) return m_Fixed.test(14);
  if (dof == disp_dwdxy) return m_Fixed.test(15);
  if (dof == fluid) return m_Fixed.test(16);
  if (dof == disp_z_1) return m_Fixed.test(17);
  if (dof == disp_z_3) return m_Fixed.test(18);
  if (dof == disp_x1_2) return m_Fixed.test(19);
  if (dof == disp_x2_2) return m_Fixed.test(20);

  return false;
}

PetscScalar cBoundaryConditionStructure::getPrescribedValue(
    eKnownDofs dof, cPoint* curPt, PetscReal omega) const {
  if (dof == disp_x1) return m_Value[0];
  if (dof == disp_x2) return m_Value[1];
  if (dof == disp_x3) return m_Value[2];
  if (dof == disp_w1) return m_Value[3];
  if (dof == disp_w2) return m_Value[4];
  if (dof == disp_w3) return m_Value[5];
  if (dof == pore0) return m_Value[6];
  if (dof == pore1) return m_Value[7];
  if (dof == pore2) return m_Value[8];
  if (dof == pore3) return m_Value[9];
  if (dof == disp_xd3) return m_Value[10];
  if (dof == disp_wd1) return m_Value[11];
  if (dof == disp_wd2) return m_Value[12];
  if (dof == disp_dwdx) return m_Value[13];
  if (dof == disp_dwdy) return m_Value[14];
  if (dof == disp_dwdxy) return m_Value[15];
  if (dof == fluid) return m_Value[16];
  if (dof == disp_z_1) return m_Value[17];
  if (dof == disp_z_3) return m_Value[18];
  if (dof == disp_x1_2) return m_Value[19];
  if (dof == disp_x2_2) return m_Value[20];

  return 0.0;
}

/*BEGIN_NO_COVERAGE*/
std::ostream& cBoundaryConditionStructure::write(std::ostream& os) const {
  os << "Nr. ";
  os.width(2);
  os << getId();
  os << "  (Struktur-Randbedingung)" << std::endl;

  for (short k = 0; k < eNumberOfDofs; k++)
    os << std::setw(3) << m_Fixed.test(k) << std::setw(5) << m_Value[k]
       << std::endl;
  return os;
}

std::istream& cBoundaryConditionStructure::read(std::istream& is) {
  short wert;

  cId::read(is);
  is >> m_Identifier;

  for (short k = 0; k < eNumberOfDofs; k++) {
    is >> wert;
    m_Fixed.set(k, wert);
    is >> m_Value[k];
  }

  return is;
}

std::ostream& cBoundaryConditionStructure::writeXml(std::ostream& os) const {
  os << "<NodeBC Type=\"structure\" Name=\"" << getIdentifier() << "\">"
     << std::endl;
  os << "  <Id>" << getId() << "</Id>" << std::endl;
  os << "  <u1>" << m_Fixed.test(0) << "</u1><valu1>" << m_Value[0]
     << "</valu1>" << std::endl;
  os << "  <u2>" << m_Fixed.test(1) << "</u2><valu2>" << m_Value[1]
     << "</valu2>" << std::endl;
  os << "  <u3>" << m_Fixed.test(2) << "</u3><valu3>" << m_Value[2]
     << "</valu3>" << std::endl;
  os << "  <w1>" << m_Fixed.test(3) << "</w1><valw1>" << m_Value[3]
     << "</valw1>" << std::endl;
  os << "  <w2>" << m_Fixed.test(4) << "</w2><valw2>" << m_Value[4]
     << "</valw2>" << std::endl;
  os << "  <w3>" << m_Fixed.test(5) << "</w3><valw3>" << m_Value[5]
     << "</valw3>" << std::endl;
  os << "  <p0>" << m_Fixed.test(6) << "</p0><valp0>" << m_Value[6]
     << "</valp0>" << std::endl;
  os << "  <p1>" << m_Fixed.test(7) << "</p1><valp1>" << m_Value[7]
     << "</valp1>" << std::endl;
  os << "  <p2>" << m_Fixed.test(8) << "</p2><valp2>" << m_Value[8]
     << "</valp2>" << std::endl;
  os << "  <p3>" << m_Fixed.test(9) << "</p3><valp3>" << m_Value[9]
     << "</valp3>" << std::endl;
  os << "  <xd3>" << m_Fixed.test(10) << "</xd3><valxd3>" << m_Value[10]
     << "</valxd3>" << std::endl;
  os << "  <wd1>" << m_Fixed.test(11) << "</wd1><valwd1>" << m_Value[11]
     << "</valwd1>" << std::endl;
  os << "  <wd2>" << m_Fixed.test(12) << "</wd2><valwd2>" << m_Value[12]
     << "</valwd2>" << std::endl;
  os << "  <dwdx>" << m_Fixed.test(13) << "</dwdx><valdwdx>" << m_Value[13]
     << "</valdwdx>" << std::endl;
  os << "  <dwdy>" << m_Fixed.test(14) << "</dwdy><valdwdy>" << m_Value[14]
     << "</valdwdy>" << std::endl;
  os << "  <dwdxy>" << m_Fixed.test(15) << "</dwdxy><valdwdxy>" << m_Value[15]
     << "</valdwdxy>" << std::endl;
  os << "  <fluid>" << m_Fixed.test(16) << "</fluid><valfluid>" << m_Value[16]
     << "</valfluid>" << std::endl;
  os << "  <disp_z_1>" << m_Fixed.test(17) << "</disp_z_1><valdisp_z_1>"
     << m_Value[17] << "</valdisp_z_1>" << std::endl;
  os << "  <disp_z_3>" << m_Fixed.test(18) << "</disp_z_3><valdisp_z_3>"
     << m_Value[18] << "</valdisp_z_3>" << std::endl;
  os << "  <disp_x1_2>" << m_Fixed.test(19) << "</disp_x1_2><valdisp_x1_2>"
     << m_Value[19] << "</valdisp_x1_2>" << std::endl;
  os << "  <disp_x2_2>" << m_Fixed.test(20) << "</disp_x2_2><valdisp_x2_2>"
     << m_Value[20] << "</valdisp_x2_2>" << std::endl;
  os << "</NodeBC>" << std::endl;

  return os;
}
/*END_NO_COVERAGE*/