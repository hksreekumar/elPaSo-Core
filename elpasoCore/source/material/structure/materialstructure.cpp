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

#include "materialstructure.h"

/*BEGIN_NO_COVERAGE*/
cMaterialStructure::cMaterialStructure()
{
  m_T  = 0.;
  m_A  = 0.;
  m_Ix = 0.;
  m_Iy = 0.;
  m_Iz = 0.;
  m_Ks = 5./6.;
  m_Fi = 0.;
}

cMaterialStructure::cMaterialStructure(const cMaterialStructure &other) :
  cMaterial(other)
{
  m_T  = other.getT();
  m_A  = other.getA();
  m_Ix = other.getIx();
  m_Iy = other.getIy();
  m_Iz = other.getIz();
  m_Ks = other.getKs();
  m_Fi = other.getFi();
}

cMaterialStructure::~cMaterialStructure()
{
  // empty
}
/*END_NO_COVERAGE*/