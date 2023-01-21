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

#include "materialspring.h"

/*BEGIN_NO_COVERAGE*/
cMaterialSpring::cMaterialSpring(void) :
  cMaterialStructureIsoLin()
{
  m_Cx = 0.;
  m_Cy = 0.;
  m_Cz = 0.;
  m_Crx = 0.;
  m_Cry = 0.;
  m_Crz = 0.;

  m_Cx_vis = 0.;
  m_Cy_vis = 0.;
  m_Cz_vis = 0.;
  m_Crx_vis = 0.;
  m_Cry_vis = 0.;
  m_Crz_vis = 0.;

  m_Eta_Spring = 0.;
  m_ViscoType = -1;
  #ifdef PETSC_USE_COMPLEX
    m_ImagOne = std::complex<PetscReal>(0., 1.);
  #else

  #endif
}


cMaterialSpring::cMaterialSpring(const cMaterialSpring &other) :
  cMaterialStructureIsoLin(other)
{
  #ifdef PETSC_USE_COMPLEX
    m_Cx = other.getCx().real();
    m_Cy = other.getCy().real();
    m_Cz = other.getCz().real();
    m_Crx = other.getCrx().real();
    m_Cry = other.getCry().real();
    m_Crz = other.getCrz().real();
    m_ImagOne = std::complex<PetscReal>(0., 1.);
  #else
    m_Cx = other.getCx();
    m_Cy = other.getCy();
    m_Cz = other.getCz();
    m_Crx = other.getCrx();
    m_Cry = other.getCry();
    m_Crz = other.getCrz();
  #endif

  m_Cx_vis = other.getCx();
  m_Cy_vis = other.getCy();
  m_Cz_vis = other.getCz();
  m_Crx_vis = other.getCrx();
  m_Cry_vis = other.getCry();
  m_Crz_vis = other.getCrz();

  m_Eta_Spring = other.getEtaSpring();
  m_ViscoType = other.getType();
}


cMaterialSpring::~ cMaterialSpring()
{
  // leer
}

void cMaterialSpring::readEtaSpringFromFile(void)
{
  std::ifstream file;
  PetscReal Frequency;

  file.open(m_Filename_Eta_Spring.c_str());

  if (file.fail())
  {
    throw cException("Error reading file containing spring values: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(file, line); ++i)
      ;
  m_lines_in_Eta_Spring_file = i;

  m_Omega_Eta_Spring.resize(m_lines_in_Eta_Spring_file);
  m_Eta_Spring_var.resize(m_lines_in_Eta_Spring_file);

  file.close();
  file.open(m_Filename_Eta_Spring.c_str());

  for (PetscInt k=0; k<m_lines_in_Eta_Spring_file; k++)
  {
    file >> Frequency;
    m_Omega_Eta_Spring.at(k) = 2. * M_PI * Frequency;
    file >> m_Eta_Spring_var.at(k);
  }

  file.close();
}

void cMaterialSpring::readCxFromFile(void)
{
  std::ifstream file;
  PetscReal Frequency;

  file.open(m_Filename_Cx.c_str());

  if (file.fail())
  {
    throw cException("Error reading file containing spring values: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(file, line); ++i)
      ;
  m_lines_in_Cx_file = i;

  m_Omega_Cx.resize(m_lines_in_Cx_file);
  m_Cx_var.resize(m_lines_in_Cx_file);

  file.close();
  file.open(m_Filename_Cx.c_str());

  for (PetscInt k=0; k<m_lines_in_Cx_file; k++)
  {
    file >> Frequency;
    m_Omega_Cx.at(k) = 2. * M_PI * Frequency;
    file >> m_Cx_var.at(k);
  }

  file.close();
}

void cMaterialSpring::readCyFromFile(void)
{
  std::ifstream file;
  PetscReal Frequency;

  file.open(m_Filename_Cy.c_str());

  if (file.fail())
  {
    throw cException("Error reading file containing spring values: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(file, line); ++i)
      ;
  m_lines_in_Cy_file = i;

  m_Omega_Cy.resize(m_lines_in_Cy_file);
  m_Cy_var.resize(m_lines_in_Cy_file);

  file.close();
  file.open(m_Filename_Cy.c_str());

  for (PetscInt k=0; k<m_lines_in_Cy_file; k++)
  {
    file >> Frequency;
    m_Omega_Cy.at(k) = 2. * M_PI * Frequency;
    file >> m_Cy_var.at(k);
  }

  file.close();
}

void cMaterialSpring::readCzFromFile(void)
{
  std::ifstream file;
  PetscReal Frequency;

  file.open(m_Filename_Cz.c_str());

  if (file.fail())
  {
    throw cException("Error reading file containing spring values: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(file, line); ++i)
      ;
  m_lines_in_Cz_file = i;

  m_Omega_Cz.resize(m_lines_in_Cz_file);
  m_Cz_var.resize(m_lines_in_Cz_file);

  file.close();
  file.open(m_Filename_Cz.c_str());

  for (PetscInt k=0; k<m_lines_in_Cz_file; k++)
  {
    file >> Frequency;
    m_Omega_Cz.at(k) = 2. * M_PI * Frequency;
    file >> m_Cz_var.at(k);
  }

  file.close();
}
void cMaterialSpring::readCrxFromFile(void)
{
  std::ifstream file;
  PetscReal Frequency;

  file.open(m_Filename_Crx.c_str());

  if (file.fail())
  {
    throw cException("Error reading file containing spring values: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(file, line); ++i)
      ;
  m_lines_in_Crx_file = i;

  m_Omega_Crx.resize(m_lines_in_Crx_file);
  m_Crx_var.resize(m_lines_in_Crx_file);

  file.close();
  file.open(m_Filename_Crx.c_str());

  for (PetscInt k=0; k<m_lines_in_Crx_file; k++)
  {
    file >> Frequency;
    m_Omega_Crx.at(k) = 2. * M_PI * Frequency;
    file >> m_Crx_var.at(k);
  }

  file.close();
}

void cMaterialSpring::readCryFromFile(void)
{
  std::ifstream file;
  PetscReal Frequency;

  file.open(m_Filename_Cry.c_str());

  if (file.fail())
  {
    throw cException("Error reading file containing spring values: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(file, line); ++i)
      ;
  m_lines_in_Cry_file = i;

  m_Omega_Cry.resize(m_lines_in_Cry_file);
  m_Cry_var.resize(m_lines_in_Cry_file);

  file.close();
  file.open(m_Filename_Cry.c_str());

  for (PetscInt k=0; k<m_lines_in_Cry_file; k++)
  {
    file >> Frequency;
    m_Omega_Cry.at(k) = 2. * M_PI * Frequency;
    file >> m_Cry_var.at(k);
  }

  file.close();
}

void cMaterialSpring::readCrzFromFile(void)
{
  std::ifstream file;
  PetscReal Frequency;

  file.open(m_Filename_Crz.c_str());

  if (file.fail())
  {
    throw cException("Error reading file containing spring values: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(file, line); ++i)
      ;
  m_lines_in_Crz_file = i;

  m_Omega_Crz.resize(m_lines_in_Crz_file);
  m_Crz_var.resize(m_lines_in_Crz_file);

  file.close();
  file.open(m_Filename_Crz.c_str());

  for (PetscInt k=0; k<m_lines_in_Crz_file; k++)
  {
    file >> Frequency;
    m_Omega_Crz.at(k) = 2. * M_PI * Frequency;
    file >> m_Crz_var.at(k);
  }

  file.close();
}

void cMaterialSpring::updateMaterial()
{
  // Update the frequency-dependent loss factor
  if (!m_Filename_Eta_Spring.empty())
  {
    if (getOmega() <= m_Omega_Eta_Spring[0]) {
      m_Eta_Spring = m_Eta_Spring_var[0];
      message("    Eta in spring material %d taken by file and set to %f.\n", getId(), m_Eta_Spring);
    }
    else if (getOmega() >= m_Omega_Eta_Spring[m_lines_in_Eta_Spring_file-1]) {
      m_Eta_Spring = m_Eta_Spring_var[m_lines_in_Eta_Spring_file-1];
      message("    Eta in spring material %d taken by file and set to %f.\n", getId(), m_Eta_Spring);
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_Eta_Spring_file; k++)
      {
        if (m_Omega_Eta_Spring[k-1] <= getOmega() && m_Omega_Eta_Spring[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Eta_Spring = ( ( (m_Eta_Spring_var[k]-m_Eta_Spring_var[k-1])*(getOmega()-m_Omega_Eta_Spring[k-1])/(m_Omega_Eta_Spring[k]-m_Omega_Eta_Spring[k-1]) ) + m_Eta_Spring_var[k-1] );
          message("    Eta in spring material %d taken by file and set to %f.\n", getId(), m_Eta_Spring);
        }
      }
    }
  }

  // Update the frequency-dependent spring values
  if (!m_Filename_Cx.empty())
  {
    if (getOmega() <= m_Omega_Cx[0]) {
      m_Cx = m_Cx_var[0];
      message("    Cx in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cx));
    }
    else if (getOmega() >= m_Omega_Cx[m_lines_in_Cx_file-1]) {
      m_Cx = m_Cx_var[m_lines_in_Cx_file-1];
      message("    Cx in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cx));
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_Cx_file; k++)
      {
        if (m_Omega_Cx[k-1] <= getOmega() && m_Omega_Cx[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Cx = ( ( (m_Cx_var[k]-m_Cx_var[k-1])*(getOmega()-m_Omega_Cx[k-1])/(m_Omega_Cx[k]-m_Omega_Cx[k-1]) ) + m_Cx_var[k-1] );
          message("    Cx in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cx));
        }
      }
    }
  }

  // Update the frequency-dependent spring values
  if (!m_Filename_Cy.empty())
  {
    if (getOmega() <= m_Omega_Cy[0]) {
      m_Cy = m_Cy_var[0];
      message("    Cy in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cy));
    }
    else if (getOmega() >= m_Omega_Cy[m_lines_in_Cy_file-1]) {
      m_Cy = m_Cy_var[m_lines_in_Cy_file-1];
      message("    Cy in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cy));
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_Cy_file; k++)
      {
        if (m_Omega_Cy[k-1] <= getOmega() && m_Omega_Cy[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Cy = ( ( (m_Cy_var[k]-m_Cy_var[k-1])*(getOmega()-m_Omega_Cy[k-1])/(m_Omega_Cy[k]-m_Omega_Cy[k-1]) ) + m_Cy_var[k-1] );
          message("    Cy in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cy));
        }
      }
    }
  }

  // Update the frequency-dependent spring values
  if (!m_Filename_Cz.empty())
  {
    if (getOmega() <= m_Omega_Cz[0]) {
      m_Cz = m_Cz_var[0];
      message("    Cz in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cz));
    }
    else if (getOmega() >= m_Omega_Cz[m_lines_in_Cz_file-1]) {
      m_Cz = m_Cz_var[m_lines_in_Cz_file-1];
      message("    Cz in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cz));
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_Cz_file; k++)
      {
        if (m_Omega_Cz[k-1] <= getOmega() && m_Omega_Cz[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Cz = ( ( (m_Cz_var[k]-m_Cz_var[k-1])*(getOmega()-m_Omega_Cz[k-1])/(m_Omega_Cz[k]-m_Omega_Cz[k-1]) ) + m_Cz_var[k-1] );
          message("    Cz in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cz));
        }
      }
    }
  }


  // Update the frequency-dependent spring values
  if (!m_Filename_Crx.empty())
  {
    if (getOmega() <= m_Omega_Crx[0]) {
      m_Crx = m_Crx_var[0];
      message("    Crx in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Crx));
    }
    else if (getOmega() >= m_Omega_Crx[m_lines_in_Crx_file-1]) {
      m_Crx = m_Crx_var[m_lines_in_Crx_file-1];
      message("    Crx in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Crx));
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_Crx_file; k++)
      {
        if (m_Omega_Crx[k-1] <= getOmega() && m_Omega_Crx[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Crx = ( ( (m_Crx_var[k]-m_Crx_var[k-1])*(getOmega()-m_Omega_Crx[k-1])/(m_Omega_Crx[k]-m_Omega_Crx[k-1]) ) + m_Crx_var[k-1] );
          message("    Crx in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Crx));
        }
      }
    }
  }

  // Update the frequency-dependent spring values
  if (!m_Filename_Cry.empty())
  {
    if (getOmega() <= m_Omega_Cry[0]) {
      m_Cry = m_Cry_var[0];
      message("    Cry in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cry));
    }
    else if (getOmega() >= m_Omega_Cry[m_lines_in_Cry_file-1]) {
      m_Cry = m_Cry_var[m_lines_in_Cry_file-1];
      message("    Cry in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cry));
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_Cry_file; k++)
      {
        if (m_Omega_Cry[k-1] <= getOmega() && m_Omega_Cry[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Cry = ( ( (m_Cry_var[k]-m_Cry_var[k-1])*(getOmega()-m_Omega_Cry[k-1])/(m_Omega_Cry[k]-m_Omega_Cry[k-1]) ) + m_Cry_var[k-1] );
          message("    Cry in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Cry));
        }
      }
    }
  }

  // Update the frequency-dependent spring values
  if (!m_Filename_Crz.empty())
  {
    if (getOmega() <= m_Omega_Crz[0]) {
      m_Crz = m_Crz_var[0];
      message("    Crz in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Crz));
    }
    else if (getOmega() >= m_Omega_Crz[m_lines_in_Crz_file-1]) {
      m_Crz = m_Crz_var[m_lines_in_Crz_file-1];
      message("    Crz in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Crz));
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_Crz_file; k++)
      {
        if (m_Omega_Crz[k-1] <= getOmega() && m_Omega_Crz[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Crz = ( ( (m_Crz_var[k]-m_Crz_var[k-1])*(getOmega()-m_Omega_Crz[k-1])/(m_Omega_Crz[k]-m_Omega_Crz[k-1]) ) + m_Crz_var[k-1] );
          message("    Crz in spring material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_Crz));
        }
      }
    }
  }

#ifdef PETSC_USE_COMPLEX
  if (m_ViscoType == 1) {
    m_Cx_vis = m_Cx * (1. + m_ImagOne * m_Eta_Spring);
    m_Cy_vis = m_Cy * (1. + m_ImagOne * m_Eta_Spring);
    m_Cz_vis = m_Cz * (1. + m_ImagOne * m_Eta_Spring);
    m_Crx_vis = m_Crx * (1. + m_ImagOne * m_Eta_Spring);
    m_Cry_vis = m_Cry * (1. + m_ImagOne * m_Eta_Spring);
    m_Crz_vis = m_Crz * (1. + m_ImagOne * m_Eta_Spring);
    std::cout << "Calculating " << getOmega()/(2*M_PI) <<" Hz ... Loss factor of spring material: "<< m_Eta_Spring <<std::endl;
    }
  else if (m_ViscoType == 2) {
    m_Cx_vis = m_Cx * (1. + m_ImagOne * getOmega() * m_Eta_Spring);
    m_Cy_vis = m_Cy * (1. + m_ImagOne * getOmega() * m_Eta_Spring);
    m_Cz_vis = m_Cz * (1. + m_ImagOne * getOmega() * m_Eta_Spring);
    m_Crx_vis = m_Crx * (1. + m_ImagOne * getOmega() * m_Eta_Spring);
    m_Cry_vis = m_Cry * (1. + m_ImagOne * getOmega() * m_Eta_Spring);
    m_Crz_vis = m_Crz * (1. + m_ImagOne * getOmega() * m_Eta_Spring);
    std::cout << "Calculating " << getOmega()/(2*M_PI) <<" Hz ... Loss factor of spring material: "<< getOmega() * m_Eta_Spring <<std::endl;
    }
  else {
    throw cException("wrong ViscoTyp - must be 1 or 2 !", __FILE__, __LINE__);
  }
#else
  m_Cx_vis = m_Cx;
  m_Cy_vis = m_Cy;
  m_Cz_vis = m_Cz;
  m_Crx_vis = m_Crx;
  m_Cry_vis = m_Cry;
  m_Crz_vis = m_Crz;
  std::cout << "Calculating " << getOmega()/(2*M_PI) <<" Hz ... Loss factor 0 due to real numbers."<< std::endl;
#endif
}

std::istream& cMaterialSpring::read(std::istream &is)
{
  cId::read(is);
  is >> m_Filename_Cx >> m_Filename_Cy >> m_Filename_Cz >> m_Filename_Crx >> m_Filename_Cry >> m_Filename_Crz >> m_Filename_Eta_Spring >> m_ViscoType; // First streamed into filename variable, later checked if it's a float or filename (char)

  // Check if Cx is given as number or file
  if (isnumber(m_Filename_Cx))
  {
    message("    Cx in spring material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Cx);
    data >> m_Cx;
    m_Filename_Cx = "";
  }
  else
  {
    message("    Cx in spring material %d interpreted as filename.\n", getId());
    readCxFromFile();
  }

  // Check if Cy is given as number or file
  if (isnumber(m_Filename_Cy))
  {
    message("    Cy in spring material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Cy);
    data >> m_Cy;
    m_Filename_Cy = "";
  }
  else
  {
    message("    Cy in spring material %d interpreted as filename.\n", getId());
    readCyFromFile();
  }

  // Check if Cz is given as number or file
  if (isnumber(m_Filename_Cz))
  {
    message("    Cz in spring material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Cz);
    data >> m_Cz;
    m_Filename_Cz = "";
  }
  else
  {
    message("    Cz in spring material %d interpreted as filename.\n", getId());
    readCzFromFile();
  }

  // Check if Crx is given as number or file
  if (isnumber(m_Filename_Crx))
  {
    message("    Crx in spring material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Crx);
    data >> m_Crx;
    m_Filename_Crx = "";
  }
  else
  {
    message("    Crx in spring material %d interpreted as filename.\n", getId());
    readCrxFromFile();
  }

  // Check if Cry is given as number or file
  if (isnumber(m_Filename_Cry))
  {
    message("    Cry in spring material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Cry);
    data >> m_Cry;
    m_Filename_Cry = "";
  }
  else
  {
    message("    Cry in spring material %d interpreted as filename.\n", getId());
    readCryFromFile();
  }

  // Check if Crz is given as number or file
  if (isnumber(m_Filename_Crz))
  {
    message("    Crz in spring material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Crz);
    data >> m_Crz;
    m_Filename_Crz = "";
  }
  else
  {
    message("    Crz in spring material %d interpreted as filename.\n", getId());
    readCrzFromFile();
  }

  // Check if Eta spring is given as number or file
  if (isnumber(m_Filename_Eta_Spring))
  {
    message("    Eta in spring material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Eta_Spring);
    data >> m_Eta_Spring;
    m_Filename_Eta_Spring = "";
  }
  else
  {
    message("    Eta in spring material %d interpreted as filename.\n", getId());
    readEtaSpringFromFile();
  }

  return is;
}


std::ostream& cMaterialSpring::write(std::ostream &os) const
{
  os << "spring material (" << getIdentifier() << ")" << std::endl;
  os << "Id..: " << getId() << std::endl;
  os << "Cx..: " << getCx() << std::endl;
  os << "Cy..: " << getCy() << std::endl;
  os << "Cz..: " << getCz() << std::endl;
  os << "Crx.: " << getCrx() << std::endl;
  os << "Cry.: " << getCry() << std::endl;
  os << "Crz.: " << getCrz() << std::endl;
  os << "Eta.: " << getEtaSpring() << std::endl;
  os << "Type.: " << getType() << std::endl;

  return os;
}


std::ostream& cMaterialSpring::writeXml(std::ostream &os) const
{
  os << "<Material Type=\"spring\" Name=\"" << getIdentifier() << "\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<Cx>" << getCx() << "</Cx>";
  os << "<Cy>" << getCy() << "</Cy>";
  os << "<Cz>" << getCz() << "</Cz>";
  os << "<Crx>" << getCrx() << "</Crx>";
  os << "<Cry>" << getCry() << "</Cry>";
  os << "<Crz>" << getCrz() << "</Crz>";
  os << "<eta>" << getEtaSpring() << "</eta>";
  os << "<type>" << getType() << "</type>";
  os << "</Material>";

  return os;
}

// Simple test if string is a float number - ref https://stackoverflow.com/questions/29169153/how-do-i-verify-a-string-is-valid-double-even-if-it-has-a-point-in-it
bool cMaterialSpring::isnumber(std::string s)
{
    int nb_point=0;
    int nb_eE_char = 0;
    int string_size = s.length();
    for (int i=0; i<string_size;i++)
    {
        if (s[i]=='.') {
            nb_point++;
        }
        else if (s[i]=='e' or s[i]=='E') {
            nb_eE_char++;
        }
        else if (!isdigit(s[i])) {
            return false;
        }
    }
    if (nb_point<=1 and nb_eE_char<=1) {
        return true;
    }
    else {
        return false;
    }
}
/*BEGIN_NO_COVERAGE*/