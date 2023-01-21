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

#ifndef ELPASO_BCFACTORY_H
#define ELPASO_BCFACTORY_H

#include <string>

class cBoundaryCondition;
class cFemParserInterface;
class ParserNodeConstraintIdsData;

//! @brief Factory for various boundary condition types
//! @author Harikrishnan Sreekumar
//! @date 06.12.2022
class cBoundaryConditionFactory {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  cBoundaryConditionFactory();
  //! @brief Destructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  ~cBoundaryConditionFactory();
  //! @brief Function to return the boundary condition
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  virtual cBoundaryCondition* createBoundaryCondition(
      std::string _bcIdentifier, int _id, cFemParserInterface* _parser);
};
#endif