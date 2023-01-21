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

#ifndef INFAM_MESH_H
#define INFAM_MESH_H

#include <algorithm>
#include <map>
#include <string>

#include "../bc/boundarycondition.h"
#include "../element/elementcoupling.h"
#include "../element/elementfem.h"
#include "../element/fluidflow/elementff.h"
#include "../element/interface/elementinterface.h"
#include "../element/load/elementload.h"
#include "../element/ncinterface/ncelementinterface.h"
#include "../material/fluid/materialfluid.h"
#include "../material/material.h"
#include "../material/structure/materialstructure.h"
#include "../misc/log/logging.h"
#include "../misc/nodalforce/nodalforce.h"
#include "../misc/nodalmoment/nodalmoment.h"
#include "../misc/nodalpressure/nodalpressure.h"
#include "../misc/nodalvalues/nodalvalues.h"
#include "group.h"
#include "node.h"

/// Forward declarations
class cConstraint;
class cConstraintMeshTie;

typedef std::map<PetscInt, cNode *> MapNodes;
typedef std::map<PetscInt, cElementFEM *> MapElements;
typedef std::map<PetscInt, cMaterial *> MapMaterials;
typedef std::map<PetscInt, cGroup *> MapGroups;
typedef std::map<PetscInt, cBoundaryCondition *> MapBoundaryConditions;
typedef std::map<PetscInt, cElementLoad *> MapElementLoads;
typedef std::map<PetscInt, cNodalForce *> MapNodalForces;
typedef std::map<PetscInt, cNodalMoment *> MapNodalMoments;
typedef std::map<PetscInt, cNodalPressure *> MapNodalPressures;
typedef std::map<PetscInt, cNodalValues *>
    MapNodalValues;  // for values from fluid flow elements
typedef std::map<PetscInt, cElementInterface *> MapInterface;
typedef std::map<PetscInt, cNCElementInterface *> MapNCInterface;
typedef std::map<PetscInt, cElementCoupling *> MapCplFemBem;
typedef std::map<PetscInt, cElementFF *> MapElementsFF;
typedef std::map<PetscInt, cConstraint *> MapConstraint;
typedef std::map<PetscInt, cConstraintMeshTie *> MapConstraintMeshTie;

typedef MapNodes::iterator ItMapNodes;
typedef MapElements::iterator ItMapElements;
typedef MapMaterials::iterator ItMapMaterials;
typedef MapGroups::iterator ItMapGroups;
typedef MapBoundaryConditions::iterator ItMapBoundaryConditions;
typedef MapElementLoads::iterator ItMapElementLoads;
typedef MapNodalForces::iterator ItMapNodalForces;
typedef MapNodalMoments::iterator ItMapNodalMoments;
typedef MapNodalPressures::iterator ItMapNodalPressures;
typedef MapNodalValues::iterator ItMapNodalValues;
typedef MapInterface::iterator ItMapInterface;
typedef MapNCInterface::iterator ItMapNCInterface;
typedef MapCplFemBem::iterator ItMapCplFemBem;
typedef MapElementsFF::iterator ItMapElementsFF;
typedef MapConstraint::iterator ItMapConstraint;

//! @brief class that holds the discretization
//! @author Dirk Clasen
//! @date 22.06.2005
class cMesh : public virtual cLogging {
 private:
  std::string m_Filename;     ///< name of the inputfile
  std::string m_Description;  ///< short problem description

  MapNodes m_MapNodes;              ///< nodes
  MapElements m_MapElements;        ///< elements
  MapElements m_MapElementsBackUp;  ///< elements backup (geoopt, testing)
  MapMaterials m_MapMaterials;      ///< materialproperties
  MapGroups m_MapGroups;            ///< groups created within MSC Patran
  MapBoundaryConditions m_MapBoundaryConditions;  ///< nodal boundaryconditions
  MapElementLoads m_MapElementLoads;              ///< elementloads
  MapNodalForces m_MapNodalForces;                ///< nodal forces
  MapNodalMoments m_MapNodalMoments;              ///< nodal moments
  MapNodalPressures m_MapNodalPressures;          ///< nodal pressures
  MapNodalValues m_MapNodalValues;  ///< nodal values (from fluid elements)
  MapInterface m_MapInterface;      ///< interface-elements (used for coupling)
  MapNCInterface m_MapNCInterface;  ///< non-conforming interface-elements (used
                                    ///< for coupling)
  MapCplFemBem m_MapCplFemBem;      ///< FEM-BEM-coupling elements
  MapElementsFF
      m_MapElementsFF;  ///< Fluid Flow elements containing data from CFD
  MapConstraint m_MapConstraint;  ///< Constraints
  MapConstraintMeshTie m_MapConstraintMeshTie;

  std::vector<PetscInt> m_VecModRedInputNodes;   ///< modred input set
  std::vector<PetscInt> m_VecModRedOutputNodes;  ///< modred output set

  std::vector<PetscInt> m_VecModRedActiveDofs;

  cGroup *ptrCurrentGroup;

 public:
  cMesh();
  ~cMesh();

  //! assigns the name of the inputfile
  void setFilename(const std::string &Filename) { m_Filename = Filename; }

  //! returns the name of the inputfile
  std::string getFilename(void) const { return m_Filename; }

  //! stores the new descriptiontext
  void setDescription(const std::string &Description) {
    m_Description = Description;
  }

  //! returns the descriptiontext
  std::string getDescription(void) const { return m_Description; }

  //! returns the largest local element id (real elements as well as interface
  //! elements)
  PetscInt getMaxElementId(void);

  //! returns the largest local group id
  PetscInt getMaxGroupId(void);

  // --------------------------------------------------------------------------
  //   NODES
  // --------------------------------------------------------------------------

  //! returns the number of nodes
  PetscInt getNumberOfNodes(void) const { return (PetscInt)m_MapNodes.size(); }

  //! inserts a new node
  void insertNode(cNode *ptrNode);

  //! delete a node
  void eraseNode(cNode *ptrNode);

  //! returns a node
  cNode *getNode(int Id);

  //! iterator pointing to first node
  ItMapNodes getFirstNode(void) { return m_MapNodes.begin(); }

  //! iterator pointing to last node
  ItMapNodes getLastNode(void) { return m_MapNodes.end(); }

  //! returns true if node in mesh - extra function as getNode exits on NULL)
  bool nodeInMesh(int Id);

  // --------------------------------------------------------------------------
  //   GROUPS
  // --------------------------------------------------------------------------

  //! insert a new group to this mesh
  void insertGroup(cGroup *ptrGroup);

  //! set a group as current
  void makeCurrentGroup(const PetscInt &Id);

  ItMapGroups getFirstGroup(void) { return m_MapGroups.begin(); }
  ItMapGroups getLastGroup(void) { return m_MapGroups.end(); }

  //! returns the pointer to group with id equal to id.
  //! If Id is not found this function returns a NULL pointer
  cGroup *getGroup(int Id);

  // --------------------------------------------------------------------------
  //   ELEMENTS
  // --------------------------------------------------------------------------

  //! returns the number of elements
  PetscInt getNumberOfElements(void) const {
    return (PetscInt)m_MapElements.size();
  }

  //! inserts a new element
  void insertElement(cElementFEM *ptrElement);

  //! delete an element
  void eraseElement(cElementFEM *ptrElement);

  //! returns an element
  cElementFEM *getElement(int Id);

  //! iterator pointing to first element
  ItMapElements getFirstElement(void) { return m_MapElements.begin(); }

  //! iterator pointing to last element
  ItMapElements getLastElement(void) { return m_MapElements.end(); }

  //! iterator pointing to ith element
  ItMapElements getElementAti(PetscInt index) {
    return m_MapElements.find(index);
  }

  /**
   * deletes some elements on the current process - used for distributed
   * computations
   */
  void purgeElements(int Purge, int Keep);

  void cropElements(std::vector<PetscInt> elem_rank);

  //! functions to backup and restore the element list
  void backupElements();
  void restoreElements();
  void eraseElementFromBackUp(int id) { m_MapElementsBackUp.erase(id); }

  // --------------------------------------------------------------------------
  //   MATERIALS
  // --------------------------------------------------------------------------

  //! number of materialsets
  PetscInt getNumberOfMaterials(void) const {
    return (PetscInt)m_MapMaterials.size();
  }

  //! insert a material
  void insertMaterial(cMaterial *ptrMaterial);

  //! returns a material
  cMaterial *getMaterial(int Id);

  //! iterator pointing to first material
  ItMapMaterials getFirstMaterial(void) { return m_MapMaterials.begin(); }

  //! iterator pointing to last material
  ItMapMaterials getLastMaterial(void) { return m_MapMaterials.end(); }

  //! Patran doesn't assign propertysets to Hex27 elements. Therefore,
  //! we have to assign a default material. This functions looks for a
  //! valid materialset in the map an returns the first appropriate
  //! materialset that is found.
  //! @return a pointer to a valid materialset
  cMaterialFluid *getDefaultMaterialFluid(void);

  //! Patran doesn't assign propertysets to Quad9 elements. Therefore,
  //! we have to assign a default material. This functions looks for a
  //! valid materialset in the map an returns the first appropriate
  //! materialset that is found.
  //! @return a pointer to a valid materialset
  cMaterialStructure *getDefaultMaterialStructure(void);

  // --------------------------------------------------------------------------
  //   NODAL BOUNDARYCONDITIONS
  // --------------------------------------------------------------------------

  //! returns the number of nodal boundaryconditions
  PetscInt getNumberOfBoundaryConditions(void) const {
    return (PetscInt)m_MapBoundaryConditions.size();
  }

  //! inserts a new boundarycondition
  void insertBoundaryCondition(cBoundaryCondition *ptrBoundaryCondition);

  //! returns a nodal boundarycondition
  cBoundaryCondition *getBoundaryCondition(int Id);

  //! iterator pointing to first boundarycondition
  ItMapBoundaryConditions getFirstBoundaryCondition(void) {
    return m_MapBoundaryConditions.begin();
  }

  //! iterator pointing to last boundarycondition
  ItMapBoundaryConditions getLastBoundaryCondition(void) {
    return m_MapBoundaryConditions.end();
  }

  // --------------------------------------------------------------------------
  //   CONSTRAINTS
  // --------------------------------------------------------------------------

  //! @brief returns the number of constraints
  //! @author Harikrishnan Sreekumar
  PetscInt getNumberOfConstraints(void) const {
    return (PetscInt)m_MapConstraint.size();
  }

  //! @brief inserts a new constraint
  //! @author Harikrishnan Sreekumar
  void insertConstraint(cConstraint *ptrConstraint);

  //! @brief returns a constraint
  //! @author Harikrishnan Sreekumar
  cConstraint *getConstraint(int Id);

  //! @brief iterator pointing to first constraint
  //! @author Harikrishnan Sreekumar
  ItMapConstraint getFirstConstraint(void) { return m_MapConstraint.begin(); }

  //! @brief iterator pointing to last constraint
  //! @author Harikrishnan Sreekumar
  ItMapConstraint getLastConstraint(void) { return m_MapConstraint.end(); }

  //! @brief inserts a new constraint
  //! @author Shreyas Guruprasad
  void insertConstraint(cConstraintMeshTie *ptrConstraintMeshTie);

  // --------------------------------------------------------------------------
  //    ELEMENTLOADS
  // --------------------------------------------------------------------------

  //! returns the number of elementloads
  PetscInt getNumberOfElementLoads(void) const {
    return (PetscInt)m_MapElementLoads.size();
  }

  //! inserts a new elementload
  void insertElementLoad(cElementLoad *ptrElementLoad);

  //! returns a single elementload
  cElementLoad *getElementLoad(int Id);

  //! iterator pointing to first elementload
  ItMapElementLoads getFirstElementLoad(void) {
    return m_MapElementLoads.begin();
  }

  //! iterator pointing to last elementload
  ItMapElementLoads getLastElementLoad(void) { return m_MapElementLoads.end(); }

  // --------------------------------------------------------------------------
  //   NODAL FORCES
  // --------------------------------------------------------------------------

  //! returns the number of nodal forces
  PetscInt getNumberOfNodalForces(void) const {
    return (PetscInt)m_MapNodalForces.size();
  }

  //! insert a nodal force
  void insertNodalForce(cNodalForce *ptrNodalForce);

  //! returns a nodal force
  cNodalForce *getNodalForce(int Id);

  //! iterator pointing to first nodal force
  ItMapNodalForces getFirstNodalForce(void) { return m_MapNodalForces.begin(); }

  //! iterator pointing to last nodal force
  ItMapNodalForces getLastNodalForce(void) { return m_MapNodalForces.end(); }

  // --------------------------------------------------------------------------
  //   NODAL MOMENTS
  // --------------------------------------------------------------------------

  //! returns the number of nodal forces
  PetscInt getNumberOfNodalMoments(void) const {
    return (PetscInt)m_MapNodalMoments.size();
  }

  //! insert a nodal moment
  void insertNodalMoment(cNodalMoment *ptrNodalMoment);

  //! returns a nodal moment
  cNodalMoment *getNodalMoment(int Id);

  //! iterator pointing to first nodal moment
  ItMapNodalMoments getFirstNodalMoment(void) {
    return m_MapNodalMoments.begin();
  }

  //! iterator pointing to last nodal moment
  ItMapNodalMoments getLastNodalMoment(void) { return m_MapNodalMoments.end(); }

  // --------------------------------------------------------------------------
  //   NODAL PRESSURES
  // --------------------------------------------------------------------------

  //! returns the number of nodal pressures
  PetscInt getNumberOfNodalPressures(void) const {
    return (PetscInt)m_MapNodalPressures.size();
  }

  //! insert a nodal pressure
  void insertNodalPressure(cNodalPressure *ptrNodalPressure);

  //! returns a nodal pressure
  cNodalPressure *getNodalPressure(int Id);

  //! iterator pointing to first nodal pressure
  ItMapNodalPressures getFirstNodalPressure(void) {
    return m_MapNodalPressures.begin();
  }

  //! iterator pointing to last nodal force
  ItMapNodalPressures getLastNodalPressure(void) {
    return m_MapNodalPressures.end();
  }

  // --------------------------------------------------------------------------
  //   NODAL VALUES (from fluid flow elements)
  // --------------------------------------------------------------------------

  //! returns the number of nodal pressures
  PetscInt getNumberOfNodalValues(void) const {
    return (PetscInt)m_MapNodalValues.size();
  }

  //! insert a nodal pressure
  void insertNodalValues(cNodalValues *ptrNodalValue);

  //! returns a nodal pressure
  cNodalValues *getNodalValues(int Id);

  //! iterator pointing to first nodal pressure
  ItMapNodalValues getFirstNodalValues(void) {
    return m_MapNodalValues.begin();
  }

  //! iterator pointing to last nodal force
  ItMapNodalValues getLastNodalValues(void) { return m_MapNodalValues.end(); }

  // -------------------------------------------------------------------------
  //   INTERFACE ELEMENTS
  // -------------------------------------------------------------------------

  //! returns the number of interfaceelements
  PetscInt getNumberOfInterfaceElements(void) const {
    return (PetscInt)m_MapInterface.size();
  }

  //! inserts an interface-element
  void insertInterfaceElement(cElementInterface *ptrElementInterface);

  //! iterator pointing to first interfaceelement
  ItMapInterface getFirstInterfaceElement(void) {
    return m_MapInterface.begin();
  }

  //! iterator pointing to last interfaceelement
  ItMapInterface getLastInterfaceElement(void) { return m_MapInterface.end(); }

  //! returns a single interfaceelement
  cElementInterface *getInterfaceElement(int Id);

  // -------------------------------------------------------------------------
  //   Non conforming INTERFACE ELEMENTS
  // -------------------------------------------------------------------------

  //! returns the number of interfaceelements
  PetscInt getNumberOfNCInterfaceElements(void) const {
    return (PetscInt)m_MapNCInterface.size();
  }

  //! inserts an interface-element
  void insertNCInterfaceElement(cNCElementInterface *ptrElementInterface);

  //! iterator pointing to first interfaceelement
  ItMapNCInterface getFirstNCInterfaceElement(void) {
    return m_MapNCInterface.begin();
  }

  //! iterator pointing to last interfaceelement
  ItMapNCInterface getLastNCInterfaceElement(void) {
    return m_MapNCInterface.end();
  }

  //! returns a single interfaceelement
  cNCElementInterface *getNCInterfaceElement(int Id);

  // -------------------------------------------------------------------------
  //   FEM BEM Coupling
  // -------------------------------------------------------------------------

  //! returns the number of interfaceelements
  PetscInt getNumberOfCplFemBemElements(void) const {
    return (PetscInt)m_MapCplFemBem.size();
  }

  //! inserts an interface-element
  void insertCplFemBemElement(cElementCoupling *ptrElement);

  //! iterator pointing to first interfaceelement
  ItMapCplFemBem getFirstCplFemBemElement(void) {
    return m_MapCplFemBem.begin();
  }

  //! iterator pointing to last interfaceelement
  ItMapCplFemBem getLastCplFemBemElement(void) { return m_MapCplFemBem.end(); }

  //! returns a single interfaceelement
  cElementCoupling *getCplFemBemElement(int Id);

  //! deletes some interfaceelements on the local process - used for
  //! parallel computations
  void purgeInterfaceElements(int Purge, int Keep);

  //! deletes some NC interfaceelements on the local process - used for
  //! parallel computations
  //! author: Harikrishnan Sreekumar
  void purgeNcInterfaceElements(int Purge, int Keep);

  void cropInterfaceElements(std::vector<PetscInt> elem_rank);

  //! deletes some non conforming interfaceelements on the local process - used
  //! for parallel computations
  void purgeNCInterfaceElements(int Purge, int Keep);

  void cropNCInterfaceElements(std::vector<PetscInt> elem_rank);

  //! deletes some couplingelements on the local process - used for
  //! parallel computations
  void purgeCouplingElements(int Purge, int Keep);

  void cropCouplingElements(std::vector<PetscInt> elem_rank);

  //! associate the interface elements to the rooms they belong.
  void associateInterfaceElementsToRoom(void);

  void checkElementSizes(void);

  // --------------------------------------------------------------------------
  //   FlUID FLOW ELEMENTS
  // --------------------------------------------------------------------------

  //! returns the number of elements
  PetscInt getNumberOfElementsFF(void) const {
    return (PetscInt)m_MapElementsFF.size();
  }

  //! inserts a new element
  void insertElementFF(cElementFF *ptrElement);

  //! delete an element
  void eraseElementFF(cElementFF *ptrElement);

  //! returns an element
  cElementFF *getElementFF(int Id);

  //! iterator pointing to first element
  ItMapElementsFF getFirstElementFF(void) { return m_MapElementsFF.begin(); }

  //! iterator pointing to last element
  ItMapElementsFF getLastElementFF(void) { return m_MapElementsFF.end(); }

  // --------------------------------------------------------------------------
  //   ModRed Data
  // --------------------------------------------------------------------------
  //! @brief Add an active dof to list
  //! @author Harikrishnan Sreekumar
  void addActiveDof(PetscInt _dofid) {
    m_VecModRedActiveDofs.push_back(_dofid);
  }

  //! @brief Add an input node to list
  //! @author Harikrishnan Sreekumar
  void addInputNode(PetscInt _nodeid) {
    m_VecModRedInputNodes.push_back(_nodeid);
  }

  //! @brief Add an output node to list
  //! @author Harikrishnan Sreekumar
  void addOutputNode(PetscInt _nodeid) {
    m_VecModRedOutputNodes.push_back(_nodeid);
  }

  //! @brief Get routine to return the active dof list
  //! @author Harikrishnan Sreekumar
  std::vector<PetscInt> &getModredActiveDofList() {
    return m_VecModRedActiveDofs;
  }

  //! @brief Get routine to return the input signal
  //! @author Harikrishnan Sreekumar
  std::vector<PetscInt> &getModredInputNodeList() {
    return m_VecModRedInputNodes;
  }

  //! @brief Get routine to return the output signal
  //! @author Harikrishnan Sreekumar
  std::vector<PetscInt> &getModredOutputNodeList() {
    return m_VecModRedOutputNodes;
  }
};

#endif
