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

#include <map>
#include <string>

class cElementFEM;

typedef void (*_func)(void*);
typedef cElementFEM* (*_funcElement)(void*);
typedef std::pair<std::string, _func> _pair;
typedef std::pair<std::string, _funcElement> _pairElement;

//! @class Handler
//! @brief Handler for HDF5 data handles to HDF5 parser routines
//! @author Harikrishnan Sreekumar
//! @date 21.10.2020

namespace HDF5 {
class Handler {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 21.10.2020
  Handler(std::string _handleName, _func _function) {
    functionpair.first = _handleName;
    functionpair.second = _function;
    free = false;
  }

  //! States if the handle is free
  bool free = true;
  //! Map of handle name to function
  _pair functionpair;
};
class ElementHandler {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 21.10.2020
  ElementHandler() {}

  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 21.10.2020
  ElementHandler(std::string _handleName, _funcElement _function) {
    functionpair.first = _handleName;
    functionpair.second = _function;
    free = false;
  }

  //! States if the handle is free
  bool free = true;
  //! Map of handle name to function
  _pairElement functionpair;
};
class MaterialHandler {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 26.04.2021
  MaterialHandler(std::string _attributename, std::string& _container,
                  std::string _dataset) {
    m_container = &_container;
    m_attributename = _attributename;
    m_dataset = _dataset;
  }
  //! @brief Method to set the value of container
  //! @param _container: value to be set
  //! @author Harikrishnan Sreekumar
  //! @date 26.04.2021
  void setContainer(std::string _container) { *m_container = _container; }
  //! Attribute name
  std::string m_attributename;
  //! Dataset name
  std::string m_dataset;

 private:
  //! Container holding the string value
  std::string* m_container;
};
class NodeLoadsHandler {
 public:
  //! @brief Constructor
  //! @author Harikrishnan Sreekumar
  //! @date 11.05.2021
  NodeLoadsHandler(std::string _attributename, double& _container,
                   std::string _dataset) {
    m_container = &_container;
    m_attributename = _attributename;
    m_dataset = _dataset;
  }
  //! @brief Method to set the value of container
  //! @param _container: value to be set
  //! @author Harikrishnan Sreekumar
  //! @date 11.05.2021
  void setContainer(double _container) { *m_container = _container; }
  //! Attribute name
  std::string m_attributename;
  //! Dataset name
  std::string m_dataset;

 private:
  //! Container holding the double value
  double* m_container;
};
}  // namespace HDF5