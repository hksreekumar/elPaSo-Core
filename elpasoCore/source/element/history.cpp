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

#include "history.h"

cHistory::cHistory()
{
  m_stress = NULL;
  m_strain = NULL;
  m_stressOld = NULL;
  m_strainOld = NULL;
  m_yielding = false;
  m_stressElastic = NULL; 
  m_strainElastic = NULL;
  m_strainPlastic = NULL;
  m_X = 0.0;
  m_epsV = 0.0;
  m_backStress = NULL;
  m_elasticValues = false;
}

cHistory::cHistory(const cHistory &other)
{
  //empty
}

cHistory::~cHistory()
{
  //empty
}


//! set stresses
void cHistory::setStress(cElementVector &stress)
{
  if(m_stress == NULL)
  {
    m_stress = new cElementVector(6);
    m_stress->setValue(0.0);
  }
  if(m_stressOld == NULL) m_stressOld = new cElementVector(6);
#ifdef PETSC_USE_COMPLEX
  trace("cHistory::setStress(cElementVector &stress) not implemented for complex analysis");
#else
  //store old stresses into stressesOld
  for(int i=0; i<6; ++i) (*m_stressOld)[i] = (*m_stress)[i];
  //store the new stresses
  for(int i=0; i<6; ++i) (*m_stress)[i]    = stress[i];  
#endif
}

//! set strains
void cHistory::setStrain(cElementVector &strain)
{
  if (strain.abs() != 0.)
  {
    if(m_strain == NULL)
    {
      m_strain = new cElementVector(6);
      m_strain->setValue(0.0);    
    }
    if(m_strainOld == NULL) 
    {
      m_strainOld = new cElementVector(6);
    }
#ifdef PETSC_USE_COMPLEX
    trace("cHistory::setStrain(cElementVector &strain) not implemented for complex analysis");
#else
    for(int i=0; i<6; ++i) (*m_strainOld)[i] = (*m_strain)[i];
    //store the new strain
    for(int i=0; i<6; ++i) (*m_strain)[i]    = strain[i];  
#endif
  }
  else
  {
    //std::cout<<"strain == 0.\n";
  }
}


//!save stress at beginning of yielding
void cHistory::storeElasticStress()
{
  if(m_stressElastic == NULL)
  {
    m_stressElastic = new cElementVector(6);
  }
#ifdef PETSC_USE_COMPLEX
  trace("cHistory::setElasticStress() not implemented for complex analysis");
#else
  //std::cout<<"Element: "<<getId()<< ": Stress set \n";
  for(int i=0; i<6; ++i) (*m_stressElastic)[i] = (*m_stress)[i];      
#endif
}

//! correct stress to get absolute stress
void cHistory::getAbsoluteStress(cElementVector &stress)
{
#ifdef PETSC_USE_COMPLEX
  trace("cHistory::getAbsoluteStress(cElementVector &stress) not implemented for complex analysis");
#else
  if (getIdealPlastic() == true)
  {
    //keep stress at constant value if behaviour is ideal plastic
    for(int i=0; i<6; ++i) stress[i] = (*m_stressElastic)[i] ;    
  }
  else 
  {
    //add stress to yield stress
    for(int i=0; i<6; ++i) stress[i] += (*m_stressElastic)[i] ;    
  }
  //std::cout<<"Absolute stress "<<(stress)[0]<<" "<<(stress)[1]<<" "<<(stress)[2]<<" "<<(stress)[3]<<" "<<(stress)[4]<<" "<<(stress)[5]<<"\n";
#endif
}

//!check for convergence
bool cHistory::checkConvergence()
{
#ifdef PETSC_USE_COMPLEX
  trace("cHistory::checkConvergence() not implemented for complex analysis");
#else
  PetscReal err = 0.00000001;

  if (m_stressOld != NULL && m_stress != NULL && (*m_stressOld).abs() != 0.)
  {    
    PetscReal dummy;
    dummy = ((*m_stress).abs()-(*m_stressOld).abs()); 
    //std::cout<<"Element: "<<getId()<<"\n";
    //std::cout<<"m_stress "<<(*m_stress).abs()<<" m_stressold "<<(*m_stressOld).abs()<<" dummy: "<<dummy<<"\n";

    if(std::fabs(dummy) < err) return true;   
  }
  else return false;
#endif
}

//!add backStress increment to overall backStress
void cHistory::setBackStress(cElementVector &stress)
{
  if(m_backStress == NULL)
  {
    m_backStress = new cElementVector(6);
    m_backStress->setValue(0.0);    
  }

  //std::cout<<"backStress before update\n"<<*m_backStress<<"\n";
  for(int i=0; i<6; ++i) (*m_backStress)[i] += stress[i] ;    
}

