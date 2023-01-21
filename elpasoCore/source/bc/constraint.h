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

#ifndef INFAM_CONSTRAINT_H
#define INFAM_CONSTRAINT_H

#include "../fedata/node.h"
#include "../misc/counter.h"
#include "../misc/id.h"

//! @brief base class for different constraints - MPC and CouplingConstraints
//! @date 28.10.2019
//! @author Harikrishnan K. Sreekumar
class cConstraint : public cId, private tCounter<cConstraint> {
 protected:
  std::string myConstraintType;  ///< type of constraint

 public:
  //! @brief Constructor
  //! @date 28.10.2019
  //! @author Harikrishnan K. Sreekumar
  cConstraint();

  //! @brief Destructor
  //! @date 28.10.2019
  //! @author Harikrishnan K. Sreekumar
  ~cConstraint();

  //! @brief Virtual Member: Add Constraint
  //! @date 29.10.2019
  //! @author Harikrishnan K. Sreekumar
  virtual void addConstraint(cNode *_masterNode, cNode *_slaveNode) = 0;

  //! @brief Finds all degrees of freedoms which are common for master and slave
  //! @param _globalMap reference to global map
  //! @date 04.11.2019
  //! @author Harikrishnan K. Sreekumar
  virtual void generateConstrainedDofs(
      std::map<PetscInt, PetscInt> &_globalMap) = 0;

  //! @brief Get the constraint type
  //! @date 29.10.2019
  //! @author Harikrishnan K. Sreekumar
  std::string getConstraintType() { return myConstraintType; }
};

#endif
