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

#ifndef ELPASO_ELEMENTLOADFACTORY_H
#define ELPASO_ELEMENTLOADFACTORY_H

#include <string>

class cElementLoad;
class ParserElementLoadsData;

//! @brief Factory for various element load types
//! @author Harikrishnan Sreekumar
//! @date 06.01.2023
class cElementLoadFactory {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.01.2023
  cElementLoadFactory();
  //! @brief Destructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.01.2023
  ~cElementLoadFactory();
  //! @brief Function to return the element load
  //! @author Harikrishnan Sreekumar
  //! @date 06.01.2023
  virtual cElementLoad* createElementLoad(ParserElementLoadsData _data);
};
#endif