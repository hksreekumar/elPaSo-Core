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

#ifndef INFAM_ELEMENT_H
#define INFAM_ELEMENT_H

#include "../fedata/node.h"
#include "../material/material.h"
#include "../misc/counter.h"
#include "../misc/gausspoints.h"
#include "../misc/id.h"
#include "../misc/log/logging.h"
#include "../shape/shapefunctions.h"
#include "./load/elementload.h"
#include "history.h"

static cShapeFunctions m_ShapeFunctions;  ///< shape functions
static cGaussPoints m_GaussPoints;        ///< Gauss points

//! @brief virtual baseclass of all types of elements
//! @author Dirk Clasen
//! @date 08.08.2005
class cElement : public cId,
                 private tCounter<cElement>,
                 public virtual cLogging {
 private:
  const int m_NumberOfNodes;   ///< number of nodes within element
  static bool m_Symmetric;     ///< use symmetric fsi formulation
  cElement* m_ElementSelfPtr;  ///< saves pointer to itself if cElement is non
                               ///< constant
  cElement*
      m_ElementParentPtr;  ///< saves pointer to parent element, e.g. when a
                           ///< disc element is created by a shell element

 protected:
  std::vector<cNode*> m_Nodes;  ///< pointer to the nodes building the element
  bool m_nl;  ///< default: false becomes true if element is linear element

 public:
  cHistory* m_history;  ///< history informations can be stored here, if needed.
                        ///< m_history is NULL per default!
                        //! @brief Constructor
  cElement(int NumberOfNodes);
  //! @brief Copy constructor
  cElement(const cElement& other);
  //! @brief Destructor
  virtual ~cElement();

  //! @brief set flag for using symmetric fsi formulation
  static void setSymmetry(bool flag) { m_Symmetric = flag; }

  //! @brief unset flag for using symmetric fsi formulation
  static bool checkSymmetry(void) { return m_Symmetric; }

  //! @brief Return number of nodes within current element
  inline int getNumberOfNodes(void) const { return m_NumberOfNodes; }

  //! @brief get the type of problem that can be solved using this element
  virtual ePhysicsType getPhysicsType(void) const = 0;

  //! @brief return cElement* pointer to this element
  inline cElement* getElementSelfPointer(void) const {
    return m_ElementSelfPtr;
  }

  //! @brief return cElement* pointer to parent element
  inline cElement* getElementParentPointer(void) const {
    return m_ElementParentPtr;
  }

  //! @brief set cElement* pointer to parent element
  void setElementParentPointer(cElement* elementPtr) {
    m_ElementParentPtr = elementPtr;
  }

  //! @brief returns element's history
  cHistory* getHistory(void) { return m_history; }

  //! @brief return a single node of element incidence
  //! @param index index of node within element's incidence  \f$ index \in
  //! [0..m_NumberOfNodes-1] \f$
  //! @return pointer to the node
  cNode* getNode(int index) const;

  //! @brief determine the maximum/minimum distance between two nodes in
  //! this element and compare them to current min/max values.
  //! Sometimes, cstGeomEps is too big - so it seems to be better to
  //! check if it fits to the mesh ...
  void checkGeomEps(PetscReal& CurMin, PetscReal& CurMax) const;

  //! @brief Inserts a new node at index into the element
  //! @param index index where to insert the node
  //! @param ptrNode points to the node
  void setNode(int index, cNode* ptrNode);

  //! @brief returns the number of instances of this object
  static size_t howMany(void) { return tCounter<cElement>::howMany(); }

  //! @brief write this object to a stream.
  //! This is only a virtual function. Its implementation is done
  //! within from cElement derived classes
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream& write(std::ostream& os) const = 0;

  //! @brief Computes the average co-ordinate of all the node of the this
  //! Element.
  //! @return cPoint with average coordinates
  cPoint getAvgofNodePoints(void);

  // bool isLinearElement()	{ return this->linearElement; }
  bool isNonlinearElement() const { return m_nl; }
  // bool isLinearMaterial()	{ return this->linearMaterial; }
  // bool isNonlinearMaterial(){ return this->nonlinearMaterial; }
};

//! @brief overloaded output operator
inline std::ostream& operator<<(std::ostream& os, const cElement& other) {
  return other.write(os);
}

#endif
