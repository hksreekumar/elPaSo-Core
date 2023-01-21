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

#ifndef INFAM_ELEMENT_COUPLING_H
#define INFAM_ELEMENT_COUPLING_H

#include "../bc/boundaryconditionstructure.h"
#include "element.h"

/* try to include SWEBEM stuff otherwise define pseudo object type */
/* in order to make the code compile */
#ifdef HAVE_SWEBEM
#include "grid.h"
#else
namespace swebem {
typedef cNode cColPoint;
}
#endif

//! @brief coupling element between FEM and SWEBEM
//! @author Dirk Clasen
//! @date 29.09.2006
class cElementCoupling : public cElement, private tCounter<cElementCoupling> {
 private:
  PetscReal m_Orientation;  ///< direction of the normal vector of the
                            ///< structural and the fluid BEM domain

  //! @brief the ids of the BEM nodes that match to the FEM nodes
  //! stored in the vector m_Nodes
  std::vector<PetscInt> m_IdsMatchingNodes;

  //! @brief BEM nodes that are attached to the FEM nodes.
  //! The location of m_Nodes[i] and m_MatchingNodes[i] is
  //! the same.
  std::vector<swebem::cColPoint *> m_MatchingBemNodes;

 public:
  //! @brief Constructor
  cElementCoupling(const short &NumberOfNodes);
  //! @brief Copy constructor
  cElementCoupling(const cElementCoupling &other);
  //! @brief Destructor
  ~cElementCoupling();

  //! @brief returns the number of instances of this objecttype
  static size_t howMany(void) { return tCounter<cElementCoupling>::howMany(); }

  //! @brief return the value of the orientation of the normal vectors
  PetscReal getOrientation(void) const { return m_Orientation; }

  //! @brief set the value for the orientation of the normal vectors
  void setOrientation(const PetscReal &ori) { m_Orientation = ori; }

  //! @brief compute the interaction between FEM and BEM
  void computeCouplingMatrix(Mat &K, Vec &F);

  //! @brief set the id of the BEM node that corresponds to the FEM node
  //! stored at m_Nodes[Index]
  void setIdMatchingBemNode(int Index, PetscInt NodeId) {
    m_IdsMatchingNodes[Index] = NodeId;
  }

  //! @brief return the id of the BEM node that corresponds to the FEM node
  //! stored at m_Nodes[Index]
  PetscInt getIdMatchingBemNode(int Index) const {
    return m_IdsMatchingNodes[Index];
  }

  //! @brief return a coupling node
  //! @param Index index in incidence table
  //! @return pointer to a node
  swebem::cColPoint *getMatchingBemNode(int Index) const;

  //! @brief insert a coupling node
  //! @param Index index in incidence table
  //! @param ptrNode pointer to a node
  void setMatchingBemNode(int Index, swebem::cColPoint *ptrNode);

  //! @brief return area of coupling element
  PetscReal getAreaOfCouplingElement(void);

  //! @brief denotes the shape of this element
  eElementShape getElementShape(void) const { return Quadrilateral; }

  //! @brief get the type of problem that can be solved using this element
  ePhysicsType getPhysicsType(void) const { return FemBem; }

  //! @brief write this object to an outputstream
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief write this object in XML format to a stream
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream &writeXml(std::ostream &os) const;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os,
                                const cElementCoupling &other) {
  return other.write(os);
}

#endif
