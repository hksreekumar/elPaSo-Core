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

#ifndef ELPASO_NODEMOMENTFACTORY_H
#define ELPASO_NODEMOMENTFACTORY_H

#include <string>

class cNodalMoment;
class ParserNodeMomentsData;

//! @brief Factory for various nodal moment types
//! @author Harikrishnan Sreekumar
//! @date 07.12.2022
class cNodalMomentFactory {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 07.12.2022
  cNodalMomentFactory();
  //! @brief Destructor
  //! @author Harikrishnan Sreekumar
  //! @date 07.12.2022
  ~cNodalMomentFactory();
  //! @brief Function to return the nodal moment
  //! @author Harikrishnan Sreekumar
  //! @date 07.12.2022
  virtual cNodalMoment* createNodalMoment(ParserNodeMomentsData _data);
};
#endif