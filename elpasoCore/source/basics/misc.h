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

#ifndef INFAM_MISC_H
#define INFAM_MISC_H

#include <string>

//! @brief additional functions
//! @author Dominik Reifer
//! @date 30.06.2018

//! @brief Simple test if string is a float number - ref
//! https://stackoverflow.com/questions/29169153/how-do-i-verify-a-string-is-valid-double-even-if-it-has-a-point-in-it
//! Here: Improved detection of signs and scientific number format
bool isnumber(std::string s);

#endif
