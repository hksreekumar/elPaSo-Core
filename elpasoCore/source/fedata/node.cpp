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

#include "node.h"
//#include "../misc/nodalpressure/nodalpressurestructure.h"
//#include "./nodalpressure/nodalpressurestructurefreq.h"

cNode::cNode(void) :
  cId(), cPoint(), cDegreesOfFreedom()
{
  ptrNodalForce = NULL;
  ptrNodalPressure1 = NULL;
  ptrNodalPressure2 = NULL;
  ptrNodalMoment = NULL;
  ptrNodalVelocity = NULL;

  setRoom( -1 );
  setGlobalSeqId( -1 );

  m_nrAdjacentElementsVector.resize(cstNumberOfStressDofs);
  for(int i=0; i<cstNumberOfStressDofs; ++i) setnrAdjacentElements( i , 0);
  
  setNumEle(0);
}

cNode::cNode(const cNode &other) :
  cId(other), cPoint(other), cDegreesOfFreedom(other)
{
  ptrNodalForce = other.getNodalForce();
  ptrNodalPressure1 = other.getNodalPressure1();
  ptrNodalPressure2 = other.getNodalPressure2();
  ptrNodalMoment = other.getNodalMoment();
  ptrNodalVelocity = other.getNodalVelocity();

  for (ItNodeBCs it = other.getFirstBC(); it != other.getLastBC(); it++)
    m_BoundaryConditions.insert( std::pair<PetscInt, cBoundaryCondition *>(it->second->getId(), it->second) );

  setRoom( other.getRoom( ) );
  setGlobalSeqId( other.getGlobalSeqId() );

  m_nrAdjacentElementsVector.resize(cstNumberOfStressDofs);
  for(int i=0; i<cstNumberOfStressDofs; ++i) setnrAdjacentElements( i , other.getnrAdjacentElements(i));
}

cNode::~cNode()
{
  ptrNodalForce = NULL;
  ptrNodalPressure1 = NULL;
  ptrNodalPressure2 = NULL;
  ptrNodalMoment = NULL;
  ptrNodalVelocity = NULL;
}


void cNode::insertBoundaryCondition(cBoundaryCondition *ptr)
{
#ifdef PETSC_USE_DEBUG
  if (ptr == NULL)
    throw cException("Zeiger NULL", __FILE__, __LINE__);
#endif

  m_BoundaryConditions.insert( std::pair<PetscInt, cBoundaryCondition *>(ptr->getId(), ptr) );
}

bool cNode::haveNodalBC(void)
{
  if ( m_BoundaryConditions.size() > 0 ) return true;
  else                                   return false;
}

void cNode::setNodalForce(cNodalForce *ptr)
{
#ifdef PETSC_USE_DEBUG
  if (ptr == NULL)
    throw cException("Zeiger NULL", __FILE__, __LINE__);
#endif

  ptrNodalForce = ptr;
}

bool cNode::haveNodalForce(void)
{
  if(ptrNodalForce!=NULL) return true;
  else                    return false;
}

void cNode::setNodalMoment(cNodalMoment *ptr)
{
#ifdef PETSC_USE_DEBUG
  if (ptr == NULL)
    throw cException("Zeiger NULL", __FILE__, __LINE__);
#endif

  ptrNodalMoment = ptr;
}

bool cNode::haveNodalMoment(void)
{
  if(ptrNodalMoment!=NULL) return true;
  else                    return false;
}

void cNode::setNodalPressure1(cNodalPressureStructure *ptr)
{
#ifdef PETSC_USE_DEBUG
  if (ptr == NULL)
    throw cException("Zeiger NULL", __FILE__, __LINE__);
#endif

  ptrNodalPressure1 = ptr;
}

void cNode::setNodalPressure2(cNodalPressureStructure *ptr)
{
#ifdef PETSC_USE_DEBUG
  if (ptr == NULL)
    throw cException("Zeiger NULL", __FILE__, __LINE__);
#endif

  ptrNodalPressure2 = ptr;
}

void cNode::setNodalVelocity(cNodalValuesVelocity *ptr)
{
#ifdef PETSC_USE_DEBUG
  if (ptr == NULL)
    throw cException("Zeiger NULL", __FILE__, __LINE__);
#endif

  ptrNodalVelocity = ptr;
}

bool cNode::haveNodalVelocity(void)
{
  if(ptrNodalVelocity!=NULL) return true;
  else                       return false;
}



std::istream& cNode::read(std::istream &is)
{
  cId::read(is);
  cPoint::read(is);

  return is;
}


std::ostream& cNode::write(std::ostream &os) const
{
  cId::write(os);
  os << "Seq. Id : " << getGlobalSeqId() << std::endl;
  cPoint::write(os);
  os << std::endl;
//   cDegreesOfFreedom::write(os);

  return os;
}


std::ostream& cNode::writeXml(std::ostream &os) const
{
  os << "<Node>";
  os << "<Id>" << getId() << "</Id>";
  os << "<x>" << (*this)[0] << "</x>";
  os << "<y>" << (*this)[1] << "</y>";
  os << "<z>" << (*this)[2] << "</z>";
  os << "</Node>";
  return os;
}
