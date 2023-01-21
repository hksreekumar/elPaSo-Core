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

#include "element.h"


// initialize static member
bool cElement::m_Symmetric = false;

//! create empty cElement object, class cFemParser creates subclass elements that also call this constructor
cElement::cElement(int NumberOfNodes) :
  m_NumberOfNodes( NumberOfNodes ),
  m_Nodes( NumberOfNodes )
{
    m_ElementSelfPtr = this;
    m_ElementParentPtr = NULL;
  for (size_t k=0; k<m_Nodes.size(); k++)
    m_Nodes[k] = NULL;

  m_nl = false;
  setSymmetry(false);
  m_history = NULL;
}

//! create constant cElement object; not used in cFemParser
cElement::cElement(const cElement &other) :
  cId(other),
  m_NumberOfNodes( other.getNumberOfNodes() ),
  m_Nodes( other.getNumberOfNodes() )
{
    //m_ElementSelfPtr = const_cast<cElement*>(&other); // Warning: unsafe
    //um_ElementParentPtr = const_cast<cElement*>(other.getElementParentPointer()); // Warning: unsafe
  for (int k=0; k<(int)m_Nodes.size(); k++)
    setNode(k, other.getNode(k));

  m_nl = other.isNonlinearElement();
  setSymmetry(other.checkSymmetry);
}


cElement::~cElement()
{
  for (int k=0; k<(int)m_Nodes.size(); k++)
    m_Nodes[k] = NULL;
}


void cElement::checkGeomEps(PetscReal &CurMin, PetscReal &CurMax) const
{
  PetscReal cur = 0.;
  for (int origin = 0; origin < getNumberOfNodes(); origin++) {
    for (int dest = origin+1; dest < getNumberOfNodes(); dest++) {
      cur = m_Nodes[origin]->distance( *m_Nodes[dest] );
      CurMax = std::max( CurMax, cur );
      CurMin = std::min( CurMin, cur );
    }
  }
}


cNode* cElement::getNode(int index) const
{
#ifdef PETSC_USE_DEBUG
  return m_Nodes.at(index);
#else
  return m_Nodes[index];
#endif
}


void cElement::setNode(int index, cNode *ptrNode)
{
#ifdef PETSC_USE_DEBUG
  if (ptrNode == NULL)
    throw cException("Pointer to node NULL!", __FILE__, __LINE__);

  m_Nodes.at(index) = ptrNode;
#else
  m_Nodes[index] = ptrNode;
#endif
}

cPoint cElement::getAvgofNodePoints(void){    
  // foreach m_Nodes get the co-ordinates
  float x=0.0,y=0.0,z=0.0;
  for (std::vector<cNode *>::iterator i = m_Nodes.begin(); i != m_Nodes.end(); ++i)
  {
    //std::cout<<"Node:"<<(**i);
	  cNode node = **i;
	  x += node[0];
    y += node[1];
    z += node[2];
  }    
  cPoint point;
  point[0] = x/m_NumberOfNodes;
  point[1] = y/m_NumberOfNodes;
  point[2] = z/m_NumberOfNodes;
  //std::cout<<"AvgPoint:"<<point;
  return point;
}
