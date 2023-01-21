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

#include "materialstructurevisco.h"
#include <string>
#include "../../../../../misc/hdf5/inputh5.h"
#include "../../../../../basics/misc.h"

/*BEGIN_NO_COVERAGE*/
cMaterialStructureVisco::cMaterialStructureVisco()
{
  E_vis = 0.;
  m_Ks = 5./6.;
  m_ImagOne = std::complex<PetscReal>(0., 1.);

  m_E   = 0.0;
  m_Eta = 0.0;
  m_Nu  = 0.0;
  m_ViscoType = -1;

  m_Filename_Eta == "";
  m_lines_in_eta_file = 0;

  m_Filename_E == "";
  m_lines_in_E_file = 0;

  m_Filename_T == "";
  m_lines_in_T_file = 0;
  m_calculatedElementTDataLength = -1;
}

cMaterialStructureVisco::cMaterialStructureVisco(const cMaterialStructureVisco &other) : 
  cMaterialStructureIsoLin(other)
{
  E_vis = other.getEVisco();
  m_Ks = 5./6.;
  m_ImagOne = std::complex<PetscReal>(0., 1.);

  m_E = other.getE();
  m_Eta = other.getEta();
  m_Nu = other.getNu();
  m_ViscoType = other.getType();

  m_Filename_Eta = other.getFilenameEta();
  m_lines_in_eta_file = 0;

  m_Filename_E = other.getFilenameE();
  m_lines_in_E_file = 0;

  m_Filename_T = other.getFilenameT();
  m_lines_in_T_file = 0;
}

cMaterialStructureVisco::~cMaterialStructureVisco()
{
  // empty
}

void cMaterialStructureVisco::readEtaFromFile(void)
{
#ifdef HAVE_HDF5
    std::vector<double>  mat_freqvariable;
    int m, n;
    cInputSingletonH5::getInstance()->readDenseDoubleMatrix(mat_freqvariable, m, n, "", "", m_Filename_Eta);

    m_lines_in_eta_file = m;
    m_Omega_Eta.resize(m);
    m_Eta_var.resize(m);

    for (PetscInt k = 0; k < m; k++)
    {
        m_Omega_Eta.at(k) = 2. * M_PI * mat_freqvariable[k * n + 0];
        m_Eta_var.at(k) = mat_freqvariable[k * n + 1];
    }
#else
  std::ifstream etaFile;
  PetscReal Frequency;
  etaFile.open(m_Filename_Eta.c_str());

  if (etaFile.fail())
  {
    throw cException("Error reading file containing loss factors: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(etaFile, line); ++i)
      ;
  m_lines_in_eta_file = i;

  m_Omega_Eta.resize(m_lines_in_eta_file);
  m_Eta_var.resize(m_lines_in_eta_file);

  etaFile.close();
  etaFile.open(m_Filename_Eta.c_str());

  for (PetscInt k=0; k<m_lines_in_eta_file; k++)
  {
    etaFile >> Frequency;
    m_Omega_Eta.at(k) = 2. * M_PI * Frequency;
    etaFile >> m_Eta_var.at(k);
  }

  etaFile.close();
#endif
}

void cMaterialStructureVisco::readEFromFile(void)
{
#ifdef HAVE_HDF5
    std::vector<double>  mat_freqvariable;
    int m, n;
    cInputSingletonH5::getInstance()->readDenseDoubleMatrix(mat_freqvariable, m, n, "", "", m_Filename_E);

    m_lines_in_E_file = m;
    m_Omega_E.resize(m);
    m_E_var.resize(m);

    for (PetscInt k = 0; k < m; k++)
    {
        m_Omega_E.at(k) = 2. * M_PI * mat_freqvariable[k * n + 0];
        m_E_var.at(k) = mat_freqvariable[k * n + 1];
    }
#else
  std::ifstream EFile;
  PetscReal Frequency;

  EFile.open(m_Filename_E.c_str());

  if (EFile.fail())
  {
    throw cException("Error reading file containing Young's moduli: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(EFile, line); ++i)
      ;
  m_lines_in_E_file = i;

  m_Omega_E.resize(m_lines_in_E_file);
  m_E_var.resize(m_lines_in_E_file);

  EFile.close();
  EFile.open(m_Filename_E.c_str());

  for (PetscInt k=0; k<m_lines_in_E_file; k++)
  {
    EFile >> Frequency;
    m_Omega_E.at(k) = 2. * M_PI * Frequency;
    EFile >> m_E_var.at(k);
  }

  EFile.close();
#endif
}

void cMaterialStructureVisco::readTFromFile(void)
{
#ifdef HAVE_HDF5
    std::vector<double>  mat_freqvariable;
    int m, n;
    cInputSingletonH5::getInstance()->readDenseDoubleMatrix(mat_freqvariable, m, n, "", "", m_Filename_T);

    m_lines_in_T_file = m;
    m_x_Tdata.resize(m);
    m_y_Tdata.resize(m);
    m_T_var.resize(m);

    for (PetscInt k = 0; k < m; k++)
    {
        m_x_Tdata.at(k) = mat_freqvariable[k * n + 0];
        m_y_Tdata.at(k) = mat_freqvariable[k * n + 1];
        m_T_var.at(k)   = mat_freqvariable[k * n + 2];
    }
#else
  std::ifstream TFile;

  TFile.open(m_Filename_T.c_str());
  if (TFile.fail())
  {
    throw cException("Error reading file containing thickness data: Unable to open file!", __FILE__, __LINE__);
  }

  std::string line;
  int i;
  for (i = 0; std::getline(TFile, line); ++i)
      ;
  m_lines_in_T_file = i;

  if (m_lines_in_T_file < 3)
  {
      throw cException("\n\nError while reading material \"visco\": At least three thickness value triples \"x y T\" are required in thickness file \n ", __FILE__, __LINE__);
  }

  m_x_Tdata.resize(m_lines_in_T_file);
  m_y_Tdata.resize(m_lines_in_T_file);
  m_T_var.resize(m_lines_in_T_file);

  TFile.close();
  TFile.open(m_Filename_T.c_str());

  for (PetscInt k=0; k<m_lines_in_T_file; k++)
  {
    TFile >> m_x_Tdata.at(k);
    TFile >> m_y_Tdata.at(k);
    TFile >> m_T_var.at(k);
  }

  TFile.close();
#endif
}

void cMaterialStructureVisco::updateMaterial()
{
#ifdef PETSC_USE_COMPLEX
  // Update the frequency-dependent loss factor
  if (!m_Filename_Eta.empty())
  {
    if (getOmega() <= m_Omega_Eta[0]) {
      m_Eta = m_Eta_var[0];
      message("    Eta in visco material %d taken by file and set to %f.\n", getId(), m_Eta);
    }
    else if (getOmega() >= m_Omega_Eta[m_lines_in_eta_file-1]) {
      m_Eta = m_Eta_var[m_lines_in_eta_file-1];
      message("    Eta in visco material %d taken by file and set to %f.\n", getId(), m_Eta);
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_eta_file; k++)
      {
        if (m_Omega_Eta[k-1] <= getOmega() && m_Omega_Eta[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_Eta = ( ( (m_Eta_var[k]-m_Eta_var[k-1])*(getOmega()-m_Omega_Eta[k-1])/(m_Omega_Eta[k]-m_Omega_Eta[k-1]) ) + m_Eta_var[k-1] );
          message("    Eta in visco material %d taken by file and set to %f.\n", getId(), m_Eta);
        }
      }
    }
  }

  // Update the frequency-dependent Young's modulus
  if (!m_Filename_E.empty())
  {
    if (getOmega() <= m_Omega_E[0]) {
      m_E = m_E_var[0];
      message("    Young's modulus in visco material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_E));
    }
    else if (getOmega() >= m_Omega_E[m_lines_in_E_file-1]) {
      m_E = m_E_var[m_lines_in_E_file-1];
      message("    Young's modulus in visco material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_E));
    }
    else
    {
      for (PetscInt k=1; k<m_lines_in_E_file; k++)
      {
        if (m_Omega_E[k-1] <= getOmega() && m_Omega_E[k] > getOmega())
        {
          // -- interpolate linear between known values
          m_E = ( ( (m_E_var[k]-m_E_var[k-1])*(getOmega()-m_Omega_E[k-1])/(m_Omega_E[k]-m_Omega_E[k-1]) ) + m_E_var[k-1] );
          message("    Young's modulus in visco material %d taken by file and set to %f.\n", getId(), PetscAbsScalar(m_E));
        }
      }
    }
  }

  //always take the first type unless otherwise specified
  if (m_ViscoType == 2) E_vis = m_E * (1. + m_ImagOne * getOmega() * m_Eta);
  else E_vis = m_E * (1. + m_ImagOne * m_Eta);
  return;

#else
  trace("  not able to evaluate viscoelastic material! ");
  trace("  recompile code including support for complex numbers");
  ExitApp();
#endif
}
/*END_NO_COVERAGE*/

void cMaterialStructureVisco::setupCb(cElementMatrix &Cb, cElement *elemPtr) const
{
#ifdef PETSC_USE_COMPLEX
  const PetscReal   h = getT(elemPtr);
  const PetscScalar B = E_vis * h*h*h / (12. * (1. - m_Nu*m_Nu));
  const PetscScalar G = E_vis / (2.*(1.+m_Nu));

  Cb(0,0) = B;
  Cb(1,1) = B;
  Cb(2,2) = G*h*h*h/12.;
  Cb(0,1) = B*m_Nu;
  Cb(1,0) = B*m_Nu;
#else
  trace("not able to evaluate viscoelastic material!");
  trace("recompile code including support for complex numbers");
  ExitApp();
#endif
}


void cMaterialStructureVisco::setupCs(cElementMatrix &Cs, cElement *elemPtr) const
{
#ifdef PETSC_USE_COMPLEX
  const PetscScalar Ghk = E_vis / (2.*(1.+m_Nu)) * getT(elemPtr) * m_Ks;

  Cs(0,0) = Ghk;
  Cs(1,1) = Ghk;
#else
  trace("not able to evaluate viscoelastic material!");
  trace("recompile code including support for complex numbers");
  ExitApp();
#endif
}


void cMaterialStructureVisco::setupCm(cElementMatrix &Cm) const
{
#ifdef PETSC_USE_COMPLEX
  const PetscScalar f = E_vis/(1.0-m_Nu*m_Nu);

  Cm(0,0) = f;
  Cm(1,1) = f;
  Cm(2,2) = f * 0.5 * (1.0-m_Nu);
  Cm(0,1) = f * m_Nu;
  Cm(1,0) = f * m_Nu;
#else
  trace("not able to evaluate viscoelastic material!");
  trace("recompile code including support for complex numbers");
  ExitApp();
#endif
}

/*BEGIN_NO_COVERAGE*/
void cMaterialStructureVisco::setupC(cElementMatrix &C) const
{
#ifdef PETSC_USE_COMPLEX
  const PetscReal   nu = getNu();
  const PetscScalar f  = E_vis/(1.0+nu)/(1.0-2.0*nu);

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
#else
  trace("not able to evaluate viscoelastic material!");
  trace("recompile code including support for complex numbers");
  ExitApp();
#endif
}
/*END_NO_COVERAGE*/

PetscReal cMaterialStructureVisco::getT(cElement *elementPtr) const
{
  if (!m_Filename_T.empty())
  {
    if (! elementPtr)
    {
       throw cException("\n\nERROR: cMaterialStructureVisco::getT() NULL element pointer: reimplementing element type necessary! Quit. \n ", __FILE__, __LINE__);
    }

    if (elementPtr->getElementParentPointer())
    {
        elementPtr = elementPtr->getElementParentPointer();
    }

    // check for previously calculated thickness data to improve performance
    // because getT() is called on every frequency step
    for(PetscInt k=0; k<m_calculatedElementTDataLength; ++k)
    {
        if (elementPtr == m_calculatedElementTData[k].first)
        {
            return m_calculatedElementTData[k].second;
        }
    }

    PetscReal elem_mean_x = calculateMeanCoordinate(elementPtr, 0);
    PetscReal elem_mean_y = calculateMeanCoordinate(elementPtr, 1);

    // array for interpolation data, stores the 3 nearest coordinates to element's mean coordinates, thickness and distance data
    // first index: T data entry; second index: x=0, y=1, T=2, r=3 (r: distance to mean element coordinates)
    PetscReal interp_dat_xyTr[3][4];
    //initialize interpolation_nodes_xyT array
    for (int k=0; k<3; ++k)
    {
        for(int i=0; i<4; ++i)
        {
            if (0==i) {interp_dat_xyTr[k][i] = m_x_Tdata.at(k);}
            if (1==i) {interp_dat_xyTr[k][i] = m_y_Tdata.at(k);}
            if (2==i) {interp_dat_xyTr[k][i] = m_T_var.at(k);}
            if (3==i) {interp_dat_xyTr[k][i] = calculateDistance(elem_mean_x, elem_mean_y, m_x_Tdata.at(k), m_y_Tdata.at(k) );}
        }
    }

    // determine three nearest points to mean element xy-coordinates
    for (int n=3; n<m_lines_in_T_file; ++n)
    {
        PetscReal currentDist= calculateDistance(elem_mean_x, elem_mean_y, m_x_Tdata.at(n), m_y_Tdata.at(n));
        int iLargestDist = 0; // determine index of largest distance to be possibly replaced by current index n
        if (interp_dat_xyTr[0][3] < interp_dat_xyTr[1][3]) { iLargestDist = 1; }
        if (interp_dat_xyTr[iLargestDist][3] < interp_dat_xyTr[2][3]) { iLargestDist = 2; }

        if (currentDist < interp_dat_xyTr[iLargestDist][3])
        {
            interp_dat_xyTr[iLargestDist][0] = m_x_Tdata.at(n);
            interp_dat_xyTr[iLargestDist][1] = m_y_Tdata.at(n);
            interp_dat_xyTr[iLargestDist][2] = m_T_var.at(n);
            interp_dat_xyTr[iLargestDist][3] = currentDist;
        }
    }

    PetscReal x = elem_mean_x;
    PetscReal y = elem_mean_y;
    PetscReal x1 = interp_dat_xyTr[0][0];
    PetscReal y1 = interp_dat_xyTr[0][1];
    PetscReal T1 = interp_dat_xyTr[0][2];
    PetscReal x2 = interp_dat_xyTr[1][0];
    PetscReal y2 = interp_dat_xyTr[1][1];
    PetscReal T2 = interp_dat_xyTr[1][2];
    PetscReal x3 = interp_dat_xyTr[2][0];
    PetscReal y3 = interp_dat_xyTr[2][1];
    PetscReal T3 = interp_dat_xyTr[2][2];

    PetscReal u1 = x2 - x1;
    PetscReal u2 = y2 - y1;
    PetscReal v1 = x3 - x1;
    PetscReal v2 = y3 - y1;
    PetscReal T;

    // case vector u and v are linearly independent
    if (u1*v2 - u2*v1 != 0.0)
    {
        PetscReal a = ( v2*(x-x1) - v1*(y-y1)) / (u1*v2 - u2*v1);
        PetscReal b = (-u2*(x-x1) - u1*(y-y1)) / (u1*v2 - u2*v1);
        T = T1 + a*(T2-T1) + b*(T3-T1);

        m_calculatedElementTData.push_back(std::make_pair(elementPtr, T));
        m_calculatedElementTDataLength = m_calculatedElementTData.size();
        //message("    T in visco material %d taken by file and set to %f.\n", getId(), T);
        return T;
    }
    // case vector u and v are linearly dependent, all three points lie on the same straight line
    else if (0.0 == u1*v2 - u2*v1)
    {
        // project point (x,y) perpendicular on the straight line and determine distances to the given points
        PetscReal unit_ux = u1 / std::sqrt(u1*u1+u2*u2);
        PetscReal unit_uy = u2 / std::sqrt(u1*u1+u2*u2);
        PetscReal s1 = std::sqrt((x-x1)*(x-x1)*unit_ux*unit_ux + (y-y1)*(y-y1)*unit_uy*unit_uy);
        PetscReal s2 = std::sqrt((x-x2)*(x-x2)*unit_ux*unit_ux + (y-y2)*(y-y2)*unit_uy*unit_uy);
        PetscReal s3 = std::sqrt((x-x3)*(x-x3)*unit_ux*unit_ux + (y-y3)*(y-y3)*unit_uy*unit_uy);
        PetscReal l1;
        PetscReal l2;
        PetscReal lT1;
        PetscReal lT2;
        PetscReal d12;

        // determine largest distance s and discard it
        if (s1 > s2 and s1 > s3)
        {
            l1 = s2;
            lT1 = T2;
            l2 = s3;
            lT2 = T3;
            d12 = std::sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
        }
        else if (s2 > s1 and s2 > s3)
        {
            l1 = s1;
            lT1 = T1;
            l2 = s3;
            lT2 = T3;
            d12 = std::sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
        }
        else
        {
            l1 = s1;
            lT1 = T1;
            l2 = s2;
            lT2 = T2;
            d12 = std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
        }

        //determine if the projected point is in between or in the extension of segment of the two interpolation points
        if (l1 <= d12 and l2 <= d12)
        {
            T = (l2*lT1 + l1*lT2) / (l1 + l2);
        }
        else if (l1 >= d12 and l1 > l2)
        {
            T = (l1*lT2 - l2*lT1) / (l1 - l2);
        }
        else if (l2 >= d12 and l2 > l1)
        {
            T = (l2*lT1 - l1*lT2) / (l2 - l1);
        }
        else if (l1 == l2)
        {
            T = (lT1 + lT2) / 2.0;
        }
        else
        {
            T = (lT1 + lT2) / 2.0;
            std::cout << "WARNING: cMaterialStructureViscoLFFrequencyTvar::getT() unexpected interpolation result, return averaged thickness (T1+T2)/2: " << T << std::endl;
        }

        m_calculatedElementTData.push_back(std::make_pair(elementPtr, T));
        m_calculatedElementTDataLength = m_calculatedElementTData.size();
        //message("    T in visco material %d taken by file and set to %f.\n", getId(), T);
        return T;
    }
    else
    {
        throw cException("\n\nERROR: cMaterialStructureVisco::getT() Thickness could not be calculated! Quit. \n ", __FILE__, __LINE__);
    }
  }
  else
  {
    return m_T;
  }
}

PetscReal cMaterialStructureVisco::calculateMeanCoordinate(cElement* elementPtr, int coord_index) const
{
   int totalnumberofnodes = elementPtr->getNumberOfNodes();
   PetscReal buildmeanvalue = 0.0;
   for (int k=0; k<totalnumberofnodes; ++k)
   {
       cNode* current_node = elementPtr->getNode(k);
       buildmeanvalue += (*current_node)[coord_index];
   }
   return (buildmeanvalue / totalnumberofnodes);
}

PetscReal cMaterialStructureVisco::calculateDistance(PetscReal x1, PetscReal y1, PetscReal x2, PetscReal y2) const
{
    PetscReal Delta_x = x2 - x1;
    PetscReal Delta_y = y2 - y1;
    PetscReal r2 = Delta_x * Delta_x + Delta_y * Delta_y;
    return std::sqrt(r2);
}

/*BEGIN_NO_COVERAGE*/
std::ostream& cMaterialStructureVisco::write(std::ostream &os) const
{
  os << "material no. " << getId() << " (visco-elastic) ("<< getIdentifier() << ")" << std::endl;
  os << "  Typ     :  " << getType() << std::endl;
  if (!m_Filename_E.empty())
  {
    os << "  E-Modul :  " << getFilenameE() << std::endl;
  }
  else
  {
    os << "  E-Modul :  " << getE() << std::endl;
  }
  if (!m_Filename_Eta.empty())
  {
    os << "  eta     :  " << getFilenameEta() << std::endl;
  }
  else
  {
    os << "  eta     :  " << getEta() << std::endl;
  }
  os << "  nu      :  " << getNu() << std::endl;
  os << "  A       :  " << getA() << std::endl;
  os << "  Ix      :  " << getIx() << std::endl;
  os << "  Iy      :  " << getIy() << std::endl;
  os << "  Iz      :  " << getIz() << std::endl;
  os << "  rho     :  " << getRho() << std::endl;
  if (!m_Filename_T.empty())
  {
    os << "  t       :  " << getFilenameT() << std::endl;
  }
  else
  {
    os << "  t       :  " << getT() << std::endl;
  }

  return os;
}


std::ostream& cMaterialStructureVisco::writeXml(std::ostream &os) const
{
  os << "<Material Type=\"visco\" Name=\"" << getIdentifier() << "\">";
  os << "<Id>" << getId() << "</Id>";
  os << "<typ>" << getType() << "</typ>";
  if (!m_Filename_E.empty())
  {
    os << "<E>" << getFilenameE() << "</E>";
  }
  else
  {
    os << "<E>" << getE() << "</E>";
  }
  if (!m_Filename_Eta.empty())
  {
    os << "<eta>" << getFilenameEta() << "</eta>";
  }
  else
  {
    os << "<eta>" << getEta() << "</eta>";
  }
  os << "<nu>" << getNu() << "</nu>";
  os << "<A>" << getA() << "</A>";
  os << "<Ix>" << getIx() << "</Ix>";
  os << "<Iy>" << getIy() << "</Iy>";
  os << "<Iz>" << getIz() << "</Iz>";
  os << "<rho>" << getRho() << "</rho>";
  if (!m_Filename_T.empty())
  {
    os << "<t>" << getFilenameT() << "</t>";
  }
  else
  {
    os << "<t>" << getT() << "</t>";
  }
  os << "</Material>" << std::endl;

  return os;
}

std::istream& cMaterialStructureVisco::read(std::istream &is)
{
  cId::read(is);

  is >> m_ViscoType;
  is >> m_Filename_E; // First streamed into filename variable, later checked if it's a float or filename (char)
  is >> m_Filename_Eta; // First streamed into filename variable, later checked if it's a float or filename (char)
  is >> m_Nu;
  is >> m_A;
  is >> m_Ix >> m_Iy >> m_Iz;
  is >> m_Rho;
  is >> m_Filename_T; // First streamed into filename variable, later checked if it's a float or filename (char)

  // Check if eta is given as number or file
  if (isnumber(m_Filename_Eta))
  {
    message("    Eta in visco material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_Eta);
    data >> m_Eta;
    m_Filename_Eta = "";
  }
  else
  {
    message("    Eta in visco material %d interpreted as filename.\n", getId());
    readEtaFromFile();
  }

  // Check if the Young's modulus is given as number or file
  if (isnumber(m_Filename_E))
  {
    message("    Young's modulus in visco material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_E);
    data >> m_E;
    m_Filename_E = "";
  }
  else
  {
    message("    Young's modulus in visco material %d interpreted as filename.\n", getId());
    readEFromFile();
  }

  // Check if the thickness is given as number or file
  if (isnumber(m_Filename_T))
  {
    message("    Thickness in visco material %d interpreted as float number.\n", getId());
    std::stringstream data(m_Filename_T);
    data >> m_T;
    m_Filename_T = "";
  }
  else
  {
    message("    Thickness in visco material %d interpreted as filename.\n", getId());
    readTFromFile();
  }

  return is;
}
/*END_NO_COVERAGE*/