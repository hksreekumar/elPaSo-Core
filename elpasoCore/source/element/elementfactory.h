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

#ifndef ELPASO_ELEMENTFACTORY_H
#define ELPASO_ELEMENTFACTORY_H

#include <string>

class cElementFEM;

//! @brief Factory for various element types
//! @author Harikrishnan Sreekumar
//! @date 06.12.2022
class cElementFactory {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  cElementFactory();
  //! @brief Destructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  ~cElementFactory();
  //! @brief Function to return the element
  //! @author Harikrishnan Sreekumar
  //! @param _elementType Element type string
  //! @date 06.12.2022
  virtual cElementFEM* createElement(std::string _elementType);
};
#endif