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

#include "misc.h"

bool isnumber(std::string s) {
  int nb_point = 0;
  int nb_eE_char = 0;
  int string_size = s.length();
  for (int i = 0; i < string_size; i++) {
    if (s[i] == '-' or s[i] == '+') {
      if (i == 0)
        continue;
      else if (s[i - 1] == 'e' or s[i - 1] == 'E')
        continue;
    }
    if (s[i] == '.')
      nb_point++;
    else if (s[i] == 'e' or s[i] == 'E')
      nb_eE_char++;
    else if (!isdigit(s[i]))
      return false;
  }
  if (nb_point <= 1 and nb_eE_char <= 1)
    return true;
  else
    return false;
}
