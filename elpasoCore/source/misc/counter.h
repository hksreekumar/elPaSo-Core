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

#ifndef INFAM_COUNTER_H
#define INFAM_COUNTER_H

#include <iostream>

//! @brief template class for counting one object's instances
//! @author Scott Meyers
//! This was published by Scott Meyers in a C++ journal. Used to check,
//! if already an object of a specific type exists (in order to perform
//! operations only once) or for setting up an element statistic
//! @date 10.04.1998
template <typename T>
class tCounter {
 private:
  static size_t count;

 public:
  tCounter() { ++count; }
  tCounter(const tCounter&) { ++count; }
  virtual ~tCounter() { --count; }

  static size_t howMany() { return count; }
};

template <typename T>
size_t tCounter<T>::count = 0;

#endif
