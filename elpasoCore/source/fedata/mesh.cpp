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

#include "mesh.h"
#include "../bc/constraint.h"

cMesh::cMesh()
{
  ptrCurrentGroup = NULL;
}


cMesh::~cMesh()
{
  // --- elements
  for (ItMapElements it = getFirstElement(); it != getLastElement(); it++)
  {
    delete it->second;
    it->second = NULL;
  }

  // --- elementloads
  for (ItMapElementLoads it = getFirstElementLoad(); it != getLastElementLoad(); it++)
  {
    delete it->second;
    it->second = NULL;
  }

  // --- materials
  for (ItMapMaterials it = getFirstMaterial(); it != getLastMaterial(); it++)
  {
    delete it->second;
    it->second = NULL;
  }

  // --- nodes
  for (ItMapNodes it = getFirstNode(); it != getLastNode(); it++)
  {
    delete it->second;
    it->second = NULL;
  }

  // --- nodalforces
  for (ItMapNodalForces it = getFirstNodalForce(); it != getLastNodalForce(); it++)
  {
    delete it->second;
    it->second = NULL;
  }

  // --- nodalforces
  for (ItMapNodalMoments it = getFirstNodalMoment(); it != getLastNodalMoment(); it++)
  {
    delete it->second;
    it->second = NULL;
  }

 // --- nodalpressures
  for (ItMapNodalPressures it = getFirstNodalPressure(); it != getLastNodalPressure(); it++)
  {
    delete it->second;
    it->second = NULL;
  }

  // --- nodal boundary conditions
  for (ItMapBoundaryConditions it = getFirstBoundaryCondition(); it != getLastBoundaryCondition(); it++)
  {
    delete it->second;
    it->second = NULL;
  }
}


PetscInt cMesh::getMaxElementId(void)
{
  PetscInt result = 0;

  if ( getNumberOfElements() > 0 )
  {
    ItMapElements  itE = getLastElement();
    itE--;
    result = itE->second->getId();
  }

  if ( getNumberOfInterfaceElements() > 0 )
  {
    ItMapInterface itI = getLastInterfaceElement(); itI--;
    result = std::max( result, itI->second->getId() );
  }

  return result;
}


PetscInt cMesh::getMaxGroupId(void)
{
  PetscInt result = 0;

  for (ItMapGroups it = getFirstGroup(); it != getLastGroup(); it++)
  {
    result = std::max( result, it->second->getId() );
  }

  return result;
}


void cMesh::insertNode(cNode *ptrNode)
{
#ifdef PETSC_USE_DEBUG
  if (ptrNode == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapNodes.insert( std::pair<PetscInt, cNode *>(ptrNode->getId(), ptrNode) );
}

void cMesh::eraseNode(cNode *ptrNode)
{
#ifdef PETSC_USE_DEBUG
  if (ptrNode == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif
  //erase element from map
  m_MapNodes.erase( ptrNode->getId() );
}

void cMesh::insertGroup(cGroup *ptrGroup)
{
#ifdef PETSC_USE_DEBUG
  if (ptrGroup == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapGroups.insert( std::pair<PetscInt, cGroup *>(ptrGroup->getId(), ptrGroup) );
}


void cMesh::insertElement(cElementFEM *ptrElement)
{
#ifdef PETSC_USE_DEBUG
  if (ptrElement == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapElements.insert( std::pair<PetscInt, cElementFEM *>(ptrElement->getId(), ptrElement) );

  // --- tell the current group that a new element was inserted
  if (ptrCurrentGroup != NULL)
    ptrCurrentGroup->insertElement( ptrElement->getId() );
}

void cMesh::eraseElement(cElementFEM *ptrElement)
{
#ifdef PETSC_USE_DEBUG
  if (ptrElement == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif
  //erase element from map
  m_MapElements.erase( ptrElement->getId() );
  //find element in group
  ItMapGroups itGroups;
  for (itGroups = getFirstGroup(); itGroups != getLastGroup(); itGroups++)
  {
    itGroups->second->eraseElement(ptrElement->getId());
  }
}

void cMesh::insertMaterial(cMaterial *ptrMaterial)
{
#ifdef PETSC_USE_DEBUG
  if (ptrMaterial == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapMaterials.insert( std::pair<PetscInt, cMaterial *>(ptrMaterial->getId(), ptrMaterial) );
}


void cMesh::insertBoundaryCondition(cBoundaryCondition *ptrBoundaryCondition)
{
#ifdef PETSC_USE_DEBUG
  if (ptrBoundaryCondition == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapBoundaryConditions.insert( std::pair<PetscInt, cBoundaryCondition *>(ptrBoundaryCondition->getId(), ptrBoundaryCondition) );
}

void cMesh::insertConstraint(cConstraint *ptrConstraint)
{
#ifdef PETSC_USE_DEBUG
  if (ptrConstraint == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapConstraint.insert( std::pair<PetscInt, cConstraint *>(ptrConstraint->getId(), ptrConstraint) );
}

void cMesh::insertConstraint(cConstraintMeshTie *ptrConstraintMeshTie)
{
#ifdef PETSC_USE_DEBUG
    if (ptrConstraintMeshTie == NULL)
        throw cException("pointer NULL", __FILE__, __LINE__);
#endif

    //m_MapConstraintMeshTie.insert(std::pair<PetscInt, cConstraintMeshTie*>(ptrConstraintMeshTie->getId(), ptrConstraintMeshTie));
}

void cMesh::insertElementLoad(cElementLoad *ptrElementLoad)
{
#ifdef PETSC_USE_DEBUG
  if (ptrElementLoad == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapElementLoads.insert( std::pair<PetscInt, cElementLoad *>(ptrElementLoad->getId(), ptrElementLoad) );
}

void cMesh::insertNodalForce(cNodalForce *ptrNodalForce)
{
#ifdef PETSC_USE_DEBUG
  if (ptrNodalForce == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapNodalForces.insert( std::pair<PetscInt, cNodalForce *>(ptrNodalForce->getId(), ptrNodalForce) );
}

void cMesh::insertNodalMoment(cNodalMoment *ptrNodalMoment)
{
#ifdef PETSC_USE_DEBUG
  if (ptrNodalMoment == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapNodalMoments.insert( std::pair<PetscInt, cNodalMoment *>(ptrNodalMoment->getId(), ptrNodalMoment) );
}

void cMesh::insertNodalPressure(cNodalPressure *ptrNodalPressure)
{
#ifdef PETSC_USE_DEBUG
  if (ptrNodalPressure == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapNodalPressures.insert( std::pair<PetscInt, cNodalPressure *>(ptrNodalPressure->getId(), ptrNodalPressure) );
}

void cMesh::insertNodalValues(cNodalValues *ptrNodalValues)
{
#ifdef PETSC_USE_DEBUG
  if (ptrNodalValues == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapNodalValues.insert( std::pair<PetscInt, cNodalValues *>(ptrNodalValues->getId(), ptrNodalValues) );
}


void cMesh::insertInterfaceElement(cElementInterface *ptrElementInterface)
{
#ifdef PETSC_USE_DEBUG
  if (ptrElementInterface == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapInterface.insert( std::pair<PetscInt, cElementInterface *>(ptrElementInterface->getId(), ptrElementInterface) );
}

void cMesh::insertNCInterfaceElement(cNCElementInterface* ptrNCElementInterface)
{
#ifdef PETSC_USE_DEBUG
    if (ptrNCElementInterface == NULL)
        throw cException("pointer NULL", __FILE__, __LINE__);
#endif

    m_MapNCInterface.insert(std::pair<PetscInt, cNCElementInterface*>(ptrNCElementInterface->getId(), ptrNCElementInterface));
}


void cMesh::insertCplFemBemElement(cElementCoupling *ptrElement)
{
#ifdef PETSC_USE_DEBUG
  if (ptrElement == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapCplFemBem.insert( std::pair<PetscInt, cElementCoupling *>(ptrElement->getId(), ptrElement) );
}


void cMesh::makeCurrentGroup(const PetscInt &Id)
{
  ItMapGroups it  = m_MapGroups.find(Id);
  ptrCurrentGroup = it->second;
}


cNode* cMesh::getNode(int Id)
{
  ItMapNodes it = m_MapNodes.find(Id);

  //if (it->second == NULL) 
  if (it == m_MapNodes.end() ) {
    message("*** no node with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

 return it->second;
}

bool cMesh::nodeInMesh(int Id)
{
  ItMapNodes it = m_MapNodes.find(Id);

  //if (it->second == NULL) 
  if (it == m_MapNodes.end() ) { return false;}
  else {return true;}
}


cElementFEM* cMesh::getElement(int Id)
{
  ItMapElements it = m_MapElements.find(Id);

  //if (it->second == NULL) 
  if (it == m_MapElements.end() ){
    // check if searched element is of type fluid flow
    ItMapElementsFF itFF = m_MapElementsFF.find(Id);
    // element also not found
    if (itFF == m_MapElementsFF.end() ) {
       message("*** no element with id %d found!\n", Id);
       throw cException("pointer NULL", __FILE__, __LINE__);
       }
  return NULL;
  }

  return it->second;
}

cGroup* cMesh::getGroup(int Id)
{
  ItMapGroups it = m_MapGroups.find(Id);

  //if (it->second == NULL)
  if (it == m_MapGroups.end() ){
    PetscPrintf(PETSC_COMM_SELF, "*** no group with id %d found!\n", Id);
    ExitApp();
  }

  return it->second;
}


cMaterial* cMesh::getMaterial(int Id)
{
  ItMapMaterials it = m_MapMaterials.find(Id);

  //if (it->second == NULL)   
  if (it == m_MapMaterials.end() ){
    message("*** no material with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}

cBoundaryCondition* cMesh::getBoundaryCondition(int Id)
{
  ItMapBoundaryConditions it = m_MapBoundaryConditions.find(Id);

  //if (it->second == NULL) 
  if (it == m_MapBoundaryConditions.end() ){
    message("*** no boundary condition with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}

cConstraint* cMesh::getConstraint(int Id)
{
  ItMapConstraint it = m_MapConstraint.find(Id);

  //if (it->second == NULL) 
  if (it == m_MapConstraint.end() ){
    message("*** no constraint with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}

cElementLoad* cMesh::getElementLoad(int Id)
{
  ItMapElementLoads it = m_MapElementLoads.find(Id);

  //if (it->second == NULL) 
  if (it == m_MapElementLoads.end() ) {
    message("*** no element load with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}

cNodalForce* cMesh::getNodalForce(int Id)
{
  ItMapNodalForces it = m_MapNodalForces.find( Id );

  //if (it->second == NULL) 
  if (it == m_MapNodalForces.end() ) {
    message("*** no nodal force with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}

cNodalMoment* cMesh::getNodalMoment(int Id)
{
  ItMapNodalMoments it = m_MapNodalMoments.find( Id );

  //if (it->second == NULL) 
  if (it == m_MapNodalMoments.end() ) {
    message("*** no nodal moment with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}

cNodalPressure* cMesh::getNodalPressure(int Id)
{
  ItMapNodalPressures it = m_MapNodalPressures.find( Id );

  //if (it->second == NULL) 
  if (it == m_MapNodalPressures.end() ) {
    message("*** no nodal pressure with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}


cElementInterface* cMesh::getInterfaceElement(int Id)
{
  ItMapInterface it = m_MapInterface.find( Id );

  //if (it->second == NULL) 
  if (it == m_MapInterface.end() ) {
    message("*** no interface element with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}

cNCElementInterface* cMesh::getNCInterfaceElement(int Id)
{
    ItMapNCInterface it = m_MapNCInterface.find(Id);

    //if (it->second == NULL) 
    if (it == m_MapNCInterface.end()) {
        message("*** no interface element with id %d found!\n", Id);
        throw cException("pointer NULL", __FILE__, __LINE__);
    }

    return it->second;
}


cElementCoupling* cMesh::getCplFemBemElement(int Id)
{
  ItMapCplFemBem it = m_MapCplFemBem.find( Id );

  //if (it->second == NULL) 
  if (it == m_MapCplFemBem.end() ) {
    message("*** no FEM-BEM interface element with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}


cMaterialFluid* cMesh::getDefaultMaterialFluid(void)
{
  cMaterialFluid *ptr = NULL;

  for (ItMapMaterials it = m_MapMaterials.begin(); it!=m_MapMaterials.end(); it++)
  {
    ptr = dynamic_cast<cMaterialFluid *>(it->second);
    if (ptr != NULL)
      return ptr;
  }

  return ptr;
}


cMaterialStructure* cMesh::getDefaultMaterialStructure(void)
{
  cMaterialStructure *ptr = NULL;

  for (ItMapMaterials it = m_MapMaterials.begin(); it!=m_MapMaterials.end(); it++)
  {
    ptr = dynamic_cast<cMaterialStructure *>(it->second);
    if (ptr != NULL)
      return ptr;
  }

  return ptr;
}


void cMesh::purgeElements(int Purge, int Keep)
{
  cElement *ptrElement;
  int size = m_MapElements.size();
  ItMapElements it;
  // --- delete elements at the beginning of the map
  for (int Counter=0; Counter<Purge; Counter++)
  {
    it=m_MapElements.begin();
    ptrElement = it->second;
    eraseElement(it->second);
    //delete ptrElement;
    ptrElement = NULL;
  }
  // --- delete elements at the end of the map
  for (int Counter=size-Purge; Counter>Keep; Counter--)
  {
    it=m_MapElements.end();
    it--;
    ptrElement=it->second;
    eraseElement(it->second);
    //delete ptrElement;
    ptrElement = NULL;
  }
}

void cMesh::cropElements(std::vector<PetscInt> elem_rank)
{
	std::sort(elem_rank.begin(), elem_rank.end());

	for (ItMapElements it = getFirstElement(); it != getLastElement(); it++) {
		cElement* ptrElement = it->second;

		if (ptrElement == NULL) {
			trace("  error: cMesh::cropElements - ptrElement == NULL ");
			return;
		}

		if (!binary_search (elem_rank.begin(), elem_rank.end(), ptrElement->getId())) {
			eraseElement(it->second);
			delete ptrElement;
			ptrElement = NULL;
		}
	}
}

void cMesh::backupElements()
{
  std::cout<<m_MapElements.size()<<" elemtens to backup\n";
  for (ItMapElements it = getFirstElement(); it != getLastElement(); it++)
  {
		cElementFEM* ptrElement = it->second;
    m_MapElementsBackUp.insert( std::pair<PetscInt, cElementFEM *>(ptrElement->getId(), ptrElement) );
  }
  std::cout<<m_MapElementsBackUp.size()<<" elemtens backuped\n";
}

void cMesh::restoreElements()
{
  m_MapElements.clear();
  std::cout<<m_MapElementsBackUp.size()<<" elemtens to restore\n";
  for (ItMapElements it = m_MapElementsBackUp.begin(); it != m_MapElementsBackUp.end(); it++)
  {
		cElementFEM* ptrElement = it->second;
    m_MapElements.insert( std::pair<PetscInt, cElementFEM *>(ptrElement->getId(), ptrElement) );
  }
  std::cout<<m_MapElements.size()<<" elemtens restored\n";
}

void cMesh::purgeInterfaceElements(int Purge, int Keep)
{
  int  Counter;
  ItMapInterface it;


  // --- delete elements at the beginning of the map
  Counter = 0;
  for (it=m_MapInterface.begin(); Counter<Purge; it++)
    Counter++;

  m_MapInterface.erase(m_MapInterface.begin(), it);

  // --- delete elements at the end of the map
  Counter = 0;
  for (it=m_MapInterface.begin(); Counter<Keep; it++)
    Counter++;
  m_MapInterface.erase(it, m_MapInterface.end());
}


void cMesh::purgeNcInterfaceElements(int Purge, int Keep)
{
    int  Counter;
    ItMapNCInterface it;


    // --- delete elements at the beginning of the map
    Counter = 0;
    for (it = m_MapNCInterface.begin(); Counter < Purge; it++)
        Counter++;

    m_MapNCInterface.erase(m_MapNCInterface.begin(), it);

    // --- delete elements at the end of the map
    Counter = 0;
    for (it = m_MapNCInterface.begin(); Counter < Keep; it++)
        Counter++;
    m_MapNCInterface.erase(it, m_MapNCInterface.end());
}

void cMesh::cropInterfaceElements(std::vector<PetscInt> elem_rank)
{
	std::sort(elem_rank.begin(), elem_rank.end());

	for (ItMapInterface it = getFirstInterfaceElement(); it != getLastInterfaceElement(); it++) {
		cElementInterface* ptrElement = it->second;

		if (ptrElement == NULL) {
			trace("  error: cMesh::cropInterfaceElements - ptrElement == NULL ");
			return;
		}

		if (!binary_search (elem_rank.begin(), elem_rank.end(), ptrElement->getId())) {
			m_MapInterface.erase(it);
			delete ptrElement;
			ptrElement = NULL;
		}
	}
}


void cMesh::purgeNCInterfaceElements(int Purge, int Keep)
{
    int  Counter;
    ItMapNCInterface it;


    // --- delete elements at the beginning of the map
    Counter = 0;
    for (it = m_MapNCInterface.begin(); Counter < Purge; it++)
        Counter++;

    m_MapNCInterface.erase(m_MapNCInterface.begin(), it);

    // --- delete elements at the end of the map
    Counter = 0;
    for (it = m_MapNCInterface.begin(); Counter < Keep; it++)
        Counter++;
    m_MapNCInterface.erase(it, m_MapNCInterface.end());
}


void cMesh::cropNCInterfaceElements(std::vector<PetscInt> elem_rank)
{
    std::sort(elem_rank.begin(), elem_rank.end());

    for (ItMapNCInterface it = getFirstNCInterfaceElement(); it != getLastNCInterfaceElement(); it++) {
        cNCElementInterface* ptrElement = it->second;

        if (ptrElement == NULL) {
            trace("  error: cMesh::cropInterfaceElements - ptrElement == NULL ");
            return;
        }

        if (!binary_search(elem_rank.begin(), elem_rank.end(), ptrElement->getId())) {
            m_MapNCInterface.erase(it);
            delete ptrElement;
            ptrElement = NULL;
        }
    }
}



void cMesh::purgeCouplingElements(int Purge, int Keep)
{
  int  Counter;
  ItMapCplFemBem it;


  // --- delete elements at the beginning of the map
  Counter = 0;
  for (it=m_MapCplFemBem.begin(); Counter<Purge; it++)
    Counter++;

  m_MapCplFemBem.erase(m_MapCplFemBem.begin(), it);

  // --- delete elements at the end of the map
  Counter = 0;
  for (it=m_MapCplFemBem.begin(); Counter<Keep; it++)
    Counter++;
  m_MapCplFemBem.erase(it, m_MapCplFemBem.end());
}


void cMesh::cropCouplingElements(std::vector<PetscInt> elem_rank)
{
	std::sort(elem_rank.begin(), elem_rank.end());

	for (ItMapCplFemBem it = getFirstCplFemBemElement(); it != getLastCplFemBemElement(); it++) {
		cElementCoupling* ptrElement = it->second;

		if (ptrElement == NULL) {
			trace("  error: cMesh::cropCouplingElements - ptrElement == NULL ");
			return;
		}

		if (!binary_search (elem_rank.begin(), elem_rank.end(), ptrElement->getId())) {
			m_MapCplFemBem.erase(it);
			delete ptrElement;
			ptrElement = NULL;
		}
	}
}


void cMesh::associateInterfaceElementsToRoom(void)
{
  PetscInt Count[3] = {0,0,0}; // 0: all  1: room A  2: room B


  ItMapInterface                      itIE;        // for iterating interface elements
  std::set<PetscInt>::const_iterator  itRooms;     // for iterating room's nodes


  trace("  associating Interface Elements to Rooms ... ");

  for (itIE = getFirstInterfaceElement(); itIE != getLastInterfaceElement(); itIE++)
  {
    Count[0]++;

    const int nnod = itIE->second->getNumberOfNodes();
    std::vector<int> room( nnod, 0 );

    for (int k=0; k<nnod; k++)
    {
      room[k] = itIE->second->getNode( k )->getRoom();
    }

    int ones   = 0;
    int twos   = 0;
    int noroom = 0;

#ifdef __sun
    std::count(room.begin(), room.end(),  1, ones);
    std::count(room.begin(), room.end(),  2, twos);
    std::count(room.begin(), room.end(), -1, noroom);
#else
    ones   =  (int)std::count(room.begin(), room.end(),  1);
    twos   =  (int)std::count(room.begin(), room.end(),  2);
    noroom =  (int)std::count(room.begin(), room.end(), -1);
#endif

    if ((ones != nnod) && (twos != nnod) && (noroom != nnod))
    {
      // error: nodes belong to no ore more than one room
      trace("  Error ! Nodes associated to no or two rooms!");
      message("    nnod   = %4d\n", nnod);
      message("    ones   = %4d\n", ones);
      message("    twos   = %4d\n", twos);
      message("    noroom = %4d\n", noroom);

      for (int k=0; k<nnod; k++)
        message("room[k] = %d ", room[k]);
      message("\n");

      MPI_Abort(PETSC_COMM_WORLD, PetscGlobalRank);
    }
    else
    {
      if (ones == nnod)
      {
        itIE->second->setRoomId( 1 );
        Count[1]++;
      }
      if (twos == nnod)
      {
        itIE->second->setRoomId( 2 );
        Count[2]++;
      }
    }
  }

  message("\n");
  message("    LOCAL Number of El.: %d\n", Count[0]);
  message("    Elements in Room A : %d\n", Count[1]);
  message("    Elements in Room B : %d\n", Count[2]);
  message("  finished!\n");
}


void cMesh::checkElementSizes(void)
{
  ItMapElements it;
  PetscReal minval = 1.0e20; // min. distance between two element nodes
  PetscReal maxval = 0.0;    // max. distance between two element nodes

  trace("checking distances between element nodes ...");

  for (it = getFirstElement(); it != getLastElement(); it++) {
    it->second->checkGeomEps(minval, maxval);
  }

  message("   cstGeomEps    : %13.5e\n", cstGeomEps);
  message("   min. distance : %13.5e\n", minval);
  message("   max. distance : %13.5e\n", maxval);

  // --- minval will be zero for spring elements
  //     therefore, check for small but not zero
  if ((minval <= cstGeomEps) && (minval != 0.0)) {
    throw cException("check your cstGeomEps in mytypes.h - it is too small!", __FILE__, __LINE__);
  }
}

void cMesh::insertElementFF(cElementFF *ptrElement)
{
#ifdef PETSC_USE_DEBUG
  if (ptrElement == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif

  m_MapElementsFF.insert( std::pair<PetscInt, cElementFF *>(ptrElement->getId(), ptrElement) );

  // --- tell the current group that a new element was inserted
  if (ptrCurrentGroup != NULL)
    ptrCurrentGroup->insertElement( ptrElement->getId() );
}

void cMesh::eraseElementFF(cElementFF *ptrElement)
{
#ifdef PETSC_USE_DEBUG
  if (ptrElement == NULL)
    throw cException("pointer NULL", __FILE__, __LINE__);
#endif
  //erase element from map
  m_MapElementsFF.erase( ptrElement->getId() );
  //find element in group
  ItMapGroups itGroups;
  for (itGroups = getFirstGroup(); itGroups != getLastGroup(); itGroups++)
  {
    itGroups->second->eraseElement(ptrElement->getId());
  }
}

cElementFF* cMesh::getElementFF(int Id)
{
  ItMapElementsFF it = m_MapElementsFF.find(Id);

  //if (it->second == NULL) 
  if (it == m_MapElementsFF.end() ){
    message("*** no element with id %d found!\n", Id);
    throw cException("pointer NULL", __FILE__, __LINE__);
  }

  return it->second;
}
