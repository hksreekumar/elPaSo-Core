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

#ifndef ELPASO_MATERIALFACTORY_H
#define ELPASO_MATERIALFACTORY_H

#include <string>

class cMaterial;
class ParserMaterialData;

//! @brief Factory for various material types
//! @author Harikrishnan Sreekumar
//! @date 06.12.2022
class cMaterialFactory {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  cMaterialFactory();
  //! @brief Destructor
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  ~cMaterialFactory();
  //! @brief Function to return the material
  //! @author Harikrishnan Sreekumar
  //! @date 06.12.2022
  virtual cMaterial* createMaterial(ParserMaterialData _matData);
};
#endif