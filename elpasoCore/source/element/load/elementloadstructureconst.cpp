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

#include "elementloadstructureconst.h"

cElementLoadStructureConst::cElementLoadStructureConst() {
  // m_F[0] = std::complex<PetscReal> (0.,0.);
  // m_F[1] = std::complex<PetscReal> (0.,0.);
  // m_F[2] = std::complex<PetscReal> (0.,0.);
  m_F[0] = 0.;
  m_F[1] = 0.;
  m_F[2] = 0.;
  m_Phase = 0.;
}

cElementLoadStructureConst::cElementLoadStructureConst(
    const cElementLoadStructureConst &other)
    : cElementLoadStructure(other) {
  m_F[0] = other.getForceComponent(0, 0);
  m_F[1] = other.getForceComponent(0, 1);
  m_F[2] = other.getForceComponent(0, 2);
}

cElementLoadStructureConst::~cElementLoadStructureConst() {
  // leer
}

std::ostream &cElementLoadStructureConst::write(std::ostream &os) const {
  cId::write(os);
  os << std::endl;
  os << " Fx = " << m_F[0] << std::endl;
  os << " Fy = " << m_F[1] << std::endl;
  os << " Fz = " << m_F[2] << std::endl;
  os << " Phase = " << m_Phase << std::endl;
  return os;
}

std::istream &cElementLoadStructureConst::read(std::istream &is) {
  cId::read(is);
  PetscReal amp_x;
  PetscReal amp_y;
  PetscReal amp_z;
  is >> amp_x >> amp_y >> amp_z >> m_Phase;
#ifdef PETSC_USE_COMPLEX
  m_F[0] = std::complex<PetscReal>(amp_x * cos(m_Phase), amp_x * sin(m_Phase));
  m_F[1] = std::complex<PetscReal>(amp_y * cos(m_Phase), amp_y * sin(m_Phase));
  m_F[2] = std::complex<PetscReal>(amp_z * cos(m_Phase), amp_z * sin(m_Phase));
#else
  m_F[0] = amp_x;
  m_F[1] = amp_y;
  m_F[2] = amp_z;
#endif
  // std::cout << "Phase " << m_Phase << std::endl;
  // std::cout << "Kraft in x: " << m_F[0] << std::endl;
  // std::cout << "Kraft in y: " << m_F[1] << std::endl;
  // std::cout << "Kraft in z: " << m_F[2] << std::endl;
  return is;
}

std::ostream &cElementLoadStructureConst::writeXml(std::ostream &os) const {
  os << "<ElemLoad Type=\"structure\">";
  os << "<Id>" << getId() << "</Id>";
#ifndef PETSC_USE_COMPLEX
  os << "<Fx>" << m_F[0] << "</Fx>";
  os << "<Fy>" << m_F[1] << "</Fy>";
  os << "<Fz>" << m_F[2] << "</Fz>";
#else
  os << "<Fx>" << m_F[0].real() << "</Fx>";
  os << "<Fy>" << m_F[1].real() << "</Fy>";
  os << "<Fz>" << m_F[2].real() << "</Fz>";
#endif
  os << "<Phase>" << m_Phase << "</Phase>";
  os << "</ElemLoad>";

  return os;
}
