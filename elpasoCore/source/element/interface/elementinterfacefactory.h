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

#ifndef ELPASO_ELEMENTINTERFACEFACTORY_H
#define ELPASO_ELEMENTINTERFACEFACTORY_H

#include <string>

class cElementInterface;
class AppInterfaceElements;

//! @brief Factory for various interface element types
//! @author Harikrishnan Sreekumar
//! @date 07.12.2022
class cElementInterfaceFactory {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 07.12.2022
  cElementInterfaceFactory();
  //! @brief Destructor
  //! @author Harikrishnan Sreekumar
  //! @date 07.12.2022
  ~cElementInterfaceFactory();
  //! @brief Function to return the interface element
  //! @author Harikrishnan Sreekumar
  //! @date 07.12.2022
  virtual bool createElementInterface(std::string Type,
                                      AppInterfaceElements _data,
                                      void* userData);

  /// @brief parse a <InterfaceKirch> element
  void sElementInterfaceKirch(AppInterfaceElements* App, void* userData);

  /// @brief parse a <InterfaceMindlin> element
  void sElementInterfaceMindlin(AppInterfaceElements* App, void* userData);

  /// @brief parse a single interface element of arbitrary type
  void parseSingleInterfaceElement(AppInterfaceElements* App,
                                   cElementInterface* ptrElement,
                                   void* userData);
};
#endif