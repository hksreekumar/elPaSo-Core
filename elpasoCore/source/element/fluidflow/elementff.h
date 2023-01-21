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

#ifndef INFAM_ELEMENTFF_H
#define INFAM_ELEMENTFF_H

#include "../element.h"

//! @brief class for elements containing fluid flow data
//! @author Silja Beck
//! @date 04.09.2011
class cElementFF : public cElement, private tCounter<cElementFF> {
 private:
  // same as element base class

 protected:
 public:
  cElementFF(int NumberOfNodes);
  cElementFF(const cElementFF &other);
  virtual ~cElementFF();

  //! @brief get the type of problem that can be solved using this element
  inline ePhysicsType getPhysicsType(void) const { return Undefined; }

  //! @brief returns the number of instances of this object
  static size_t howMany(void) { return tCounter<cElementFF>::howMany(); }

  //! @brief write this object to a stream.
  //! This is only a virtual function. Its implementation is done
  //! within from cElement derived classes
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream &write(std::ostream &os) const;

  //! @brief write this object in XML
  //! @param os outputstream
  //! @return modified outputstream
  virtual std::ostream &writeXml(std::ostream &os) const = 0;
};

#endif
