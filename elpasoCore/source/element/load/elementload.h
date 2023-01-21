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

#ifndef INFAM_ELEMENT_LOAD_H
#define INFAM_ELEMENT_LOAD_H

#include "../../misc/counter.h"
#include "../../misc/id.h"
#include "../../misc/log/logging.h"

//! @brief virtual base class for element loads
//! @author Dirk Clasen
//! @date 31.05.2005
class cElementLoad : public cId,
                     private tCounter<cElementLoad>,
                     public virtual cLogging {
 public:
  //! @brief Constructor
  cElementLoad();
  //! @brief Copy constructor
  cElementLoad(const cElementLoad &other);
  //! @brief Destructor
  virtual ~cElementLoad();

  //! @brief return the number of objects of this type
  static size_t howMany(void) { return tCounter<cElementLoad>::howMany(); }

  //! @brief write this object to a stream
  //! @param os output stream
  //! @return modified outputstream
  virtual std::ostream &write(std::ostream &os) const = 0;

  //! @brief read data to this object from a stream
  //! @param is input stream
  //! @return modified inputstream
  virtual std::istream &read(std::istream &is) = 0;

  //! @brief write this object to XML stream
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

//! @brief overloaded outputoperator
inline std::ostream &operator<<(std::ostream &os, const cElementLoad &other) {
  return other.write(os);
}

//! @brief overloaded inputoperator
inline std::istream &operator>>(std::istream &is, cElementLoad &other) {
  return other.read(is);
}

#endif
