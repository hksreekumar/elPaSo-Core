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

#include "materialstructureisotrop.h"

/*BEGIN_NO_COVERAGE*/
cMaterialStructureIsotrop::cMaterialStructureIsotrop(void) :
  cMaterialStructureIsoLin()
{
  m_E  = 0.;
  m_Nu = 0.;
}


cMaterialStructureIsotrop::cMaterialStructureIsotrop(const cMaterialStructureIsotrop &other) :
  cMaterialStructureIsoLin(other)
{
  m_E  = other.getE();
  m_Nu = other.getNu();
}


cMaterialStructureIsotrop::~ cMaterialStructureIsotrop()
{
  // leer
}


std::istream& cMaterialStructureIsotrop::read(std::istream &is)
{
  cId::read(is);
  is >> m_E >> m_Nu >> m_A >> m_Ix >> m_Iy >> m_Iz >> m_Rho >> m_T >> m_Fi;

  return is;
}


std::ostream& cMaterialStructureIsotrop::write(std::ostream &os) const
{
  os << "isotropic material (" << getIdentifier() << ")" << std::endl;
  os << "Id..: " << getId() << std::endl;
  os << "E...: " << getE() << std::endl;
  os << "ks..: " << getKs() << std::endl;
  os << "nu..: " << getNu() << std::endl;
  os << "A ..: " << getA() << std::endl;
  os << "Ix .: " << getIx() << std::endl;
  os << "Iy .: " << getIy() << std::endl;
  os << "Iz .: " << getIz() << std::endl;
  os << "rho.: " << getRho() << std::endl;
  os << "t...: " << getT() << std::endl;
  os << "Fi..: " << getFi() << std::endl;

  return os;
}


std::ostream& cMaterialStructureIsotrop::writeXml(std::ostream &os) const
{
  os << "<Material Type=\"isotropic\" Name=\"" << getIdentifier() << "\">";
  os << "<Id>" << getId() << "</Id>";
#ifdef PETSC_USE_COMPLEX
  os << "<E>" << getE().real() << "</E>";
#else
  os << "<E>" << getE() << "</E>";
#endif
  os << "<nu>" << getNu() << "</nu>";
  os << "<A>" << getA() << "</A>";
  os << "<Ix>" << getIx() << "</Ix>";
  os << "<Iy>" << getIy() << "</Iy>";
  os << "<Iz>" << getIz() << "</Iz>";
  os << "<rho>" << getRho() << "</rho>";
  os << "<t>" << getT() << "</t>";
  os << "<Fi>" << getFi() << "</Fi>";
  os << "</Material>";

  return os;
}
/*END_NO_COVERAGE*/

void cMaterialStructureIsotrop::setupCb(cElementMatrix &Cb, cElement *elemPtr) const
{
  const PetscReal   h = getT();
  const PetscScalar B = m_E * h*h*h / (12. * (1. - m_Nu*m_Nu));
  const PetscScalar G = m_E / (2.*(1.+m_Nu));

  Cb(0,0) = B;
  Cb(1,1) = B;
  Cb(2,2) = G*h*h*h/12.;
  Cb(0,1) = B*m_Nu;
  Cb(1,0) = B*m_Nu;
}


void cMaterialStructureIsotrop::setupCs(cElementMatrix &Cs, cElement *elemPtr) const
{
  const PetscScalar Ghk = m_E / (2.*(1.+m_Nu)) * getT() * m_Ks;
  Cs(0,0) = Ghk;
  Cs(1,1) = Ghk;
}


void cMaterialStructureIsotrop::setupCm(cElementMatrix &Cm) const
{
  const PetscScalar f = m_E/(1.0-m_Nu*m_Nu);

  Cm(0,0) = f;
  Cm(1,1) = f;
  Cm(2,2) = f * 0.5 * (1.0-m_Nu);
  Cm(0,1) = f * m_Nu;
  Cm(1,0) = f * m_Nu;
}

/*BEGIN_NO_COVERAGE*/
void cMaterialStructureIsotrop::setupC(cElementMatrix &C) const
{
  const PetscReal   nu = getNu();
  const PetscScalar f  = getE()/(1.0+nu)/(1.0-2.0*nu);

  C(0,0) = f*(1.-nu);
  C(1,1) = f*(1.-nu);
  C(2,2) = f*(1.-nu);
  C(3,3) = f*(1.0-2.0*nu)/2.0;
  C(4,4) = f*(1.0-2.0*nu)/2.0;
  C(5,5) = f*(1.0-2.0*nu)/2.0;
  C(0,1) = f*nu;
  C(1,0) = f*nu;
  C(0,2) = f*nu;
  C(2,0) = f*nu;
  C(1,2) = f*nu;
  C(2,1) = f*nu;
}

void cMaterialStructureIsotrop::setupCps(cElementMatrix &Cps) const
{
  const PetscReal   nu = getNu();
  const PetscScalar f  = getE()/(1.0-nu*nu);

  Cps(0,0) = f;
  Cps(1,1) = f;
  Cps(0,1) = f*nu;
  Cps(1,0) = f*nu;
  Cps(2,2) = f*(1.0-nu)/2.;
}
/*END_NO_COVERAGE*/