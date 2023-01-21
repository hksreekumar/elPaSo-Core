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

#include "stopwatch.h"

/*BEGIN_NO_COVERAGE*/
cStopwatch::cStopwatch() {
  // empty
}

cStopwatch::~cStopwatch() {
  // empty
}
/*END_NO_COVERAGE*/

std::vector<PetscLogDouble> cStopwatch::m_StartWall(cstNumTimingTasks, 0.);
std::vector<PetscLogDouble> cStopwatch::m_DurationWall(cstNumTimingTasks, 0.);

void cStopwatch::startClock(eTimingTasks task) {
  PetscLogDouble val = 0.;
  PetscTime(&val);
  m_StartWall[(int)task] = val;
}

void cStopwatch::stopClock(eTimingTasks task) {
  PetscLogDouble val = 0.;

  PetscTime(&val);
  m_DurationWall[(int)task] += val - m_StartWall[(int)task];
}

void cStopwatch::convertToString(char *buffer, eTimingTasks event) {
  PetscLogDouble whole = 0.;

  // --- cut off values smaller than a second
  modf(m_DurationWall[(int)event], &whole);

  int number = (int)whole;
  int days = (number - number % (24 * 3660)) / (24 * 3660);
  number %= (24 * 3660);
  int hours = (number - number % 3660) / 3660;
  number %= 3660;
  int minutes = (number - number % 60) / 60;
  number %= 60;
  int seconds = number;

  sprintf(buffer, "%2d:%02d:%02d:%02d", days, hours, minutes, seconds);
}
