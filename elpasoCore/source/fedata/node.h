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

#ifndef INFAM_NODE_H
#define INFAM_NODE_H

#include "../bc/boundarycondition.h"
#include "../misc/counter.h"
#include "../misc/id.h"
#include "../misc/nodalforce/nodalforce.h"
#include "../misc/nodalmoment/nodalmoment.h"
#include "../misc/nodalpressure/nodalpressurestructure.h"
#include "../misc/nodalvalues/nodalvaluesvelocity.h"
#include "../misc/point.h"
#include "dof.h"

typedef std::map<PetscInt, cBoundaryCondition *>::const_iterator ItNodeBCs;

//! @brief nodes for finite elements
//!
//! @author Dirk Clasen, add ons by Silja Beck and Marco Schauer
//! @date 22.06.2005
class cNode : public cId,
              public cPoint,
              public cDegreesOfFreedom,
              private tCounter<cNode> {
 private:
  PetscInt m_GlobalSeqNumber;
  cNodalForce *ptrNodalForce;    ///< pointer to nodal force  (read from parser)
  cNodalMoment *ptrNodalMoment;  ///< pointer to nodal moment (read from parser)
  cNodalPressureStructure *ptrNodalPressure1;  ///< pointer to single nodal
                                               ///< pressure (read from parser)
  cNodalPressureStructure
      *ptrNodalPressure2;  ///< pointer to single nodal pressure (read from
                           ///< NodalPressureStructureFreq)
  cNodalValuesVelocity *ptrNodalVelocity;  ///< pointer to velocities from fluid
                                           ///< flow (read from NodalValues)

  std::vector<PetscInt>
      m_nrAdjacentElementsVector;  ///< nr of adjacent elements for the
                                   ///< eStresses

  short m_Room;    ///< Id of the room this node belongs to (-1 if it belongs to
                   ///< no room)
  short m_NumEle;  ///< Number of elements attached to this node

  std::map<PetscInt, cBoundaryCondition *>
      m_BoundaryConditions;  ///< boundary conditions assigned to node

 public:
  cNode(void);
  cNode(const cNode &other);
  ~cNode();

  //! all nodes can have arbitrary id i.e. they are not numbered from 1 to nnod.
  //! This numbering will be needed for stress computation at the nodes.
  inline PetscInt getGlobalSeqId(void) const { return m_GlobalSeqNumber; }

  //! assign a new node id used to address vector entries when computing
  //! stresses.
  void setGlobalSeqId(const PetscInt &number) { m_GlobalSeqNumber = number; }

  //! important for the stress computation at the nodes
  inline PetscInt getnrAdjacentElements(const PetscInt eStresses_nr) const {
    return m_nrAdjacentElementsVector[eStresses_nr];
  }

  inline void setnrAdjacentElements(const PetscInt eStresses_nr,
                                    const PetscInt number) {
    m_nrAdjacentElementsVector[eStresses_nr] = number;
  }

  inline void addAdjacentElement(const PetscInt eStresses_nr) {
    m_nrAdjacentElementsVector[eStresses_nr]++;
  }

  //! return iterator to first boundary condition assigned to this node
  inline ItNodeBCs getFirstBC(void) const {
    return m_BoundaryConditions.begin();
  }

  //! return iterator to last boundary condition assigned to this node
  inline ItNodeBCs getLastBC(void) const { return m_BoundaryConditions.end(); }

  //! return the number of boundary conditions assigned to this node
  inline PetscInt getNumberOfBCs(void) const {
    return (PetscInt)m_BoundaryConditions.size();
  }

  //! set id of the room this node is associated to
  inline void setRoom(short room) { m_Room = room; }

  //! get id of the room this node is associated to
  short getRoom(void) const { return m_Room; }

  //! return the number of instances of this object type
  static size_t howMany(void) { return tCounter<cNode>::howMany(); }

  //! assign a boundary condition to this node
  //! @param ptr pointer to boundary condition
  void insertBoundaryCondition(cBoundaryCondition *ptr);

  //! check whether node have nodal bc
  //! returns true if nodal force is set and false if not
  bool haveNodalBC(void);

  //! assign a nodal force to this node
  //! @param ptr pointer to nodal force
  void setNodalForce(cNodalForce *ptr);

  //! read nodal force applied to this node
  //! @return pointer to nodal force (NULL, if no force applied to this node)
  inline cNodalForce *getNodalForce(void) const { return ptrNodalForce; }

  //! check whether node have nodal force
  //! returns true if nodal force is set and false if not
  bool haveNodalForce(void);

  //! assign a nodal moment to this node
  //! @param ptr pointer to nodal moment
  void setNodalMoment(cNodalMoment *ptr);

  //! read nodal moment applied to this node
  //! @return pointer to nodal moment (NULL, if no force applied to this node)
  inline cNodalMoment *getNodalMoment(void) const { return ptrNodalMoment; }

  //! check whether node have nodal moment
  //! returns true if nodal moment is set and false if not
  bool haveNodalMoment(void);

  //! assign a nodal pressure to this node
  //! @param ptr pointer to nodal pressure
  void setNodalPressure1(cNodalPressureStructure *ptr);
  void setNodalPressure2(cNodalPressureStructure *ptr);

  //! read nodal pressure applied to this node
  //! @return pointer to nodal pressure (NULL, if no pressure applied to this
  //! node)
  inline cNodalPressureStructure *getNodalPressure1(void) const {
    return ptrNodalPressure1;
  }
  inline cNodalPressureStructure *getNodalPressure2(void) const {
    return ptrNodalPressure2;
  }

  //! assign a nodal velocities (from fluid flow) to this node
  //! @param ptr pointer to nodal velocities
  void setNodalVelocity(cNodalValuesVelocity *ptr);

  //! read nodal velocities applied to this node
  //! @return pointer to nodal velocities (NULL, if no velocities applied to
  //! this node)
  inline cNodalValuesVelocity *getNodalVelocity(void) const {
    return ptrNodalVelocity;
  }

  //! check whether node has nodal velocities assigned
  //! returns true if nodal velocity is set and false if not
  bool haveNodalVelocity(void);

  //! set number of elements attached to the node
  void setNumEle(int val) { m_NumEle = val; }
  //! returns number of elements attached to the node
  short getNumEle() { return m_NumEle; }
  //! counts number of elements attached to this element
  void countElement() { m_NumEle += 1; }

  //! reads a single object of a stream
  //! @param is inputstream
  //! @return modified outputstream
  std::istream &read(std::istream &is);

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cNode &other) {
  return other.read(is);
}

//! overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cNode &other) {
  return other.write(os);
}

#endif
