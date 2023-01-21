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

#ifndef INFAM_STOPWATCH_H
#define INFAM_STOPWATCH_H

#include <iostream>
#include <vector>

#include "petsctime.h"

//! The total number of timing tasks
const int cstNumTimingTasks = 1000;

//! @brief all events used for measure time used by code
enum eTimingTasks {
  Overall = 0,
  SetupKM = 1,
  Solver = 2,
  reorder = 3,
  output = 4,
  modred_overall = 11,
  modred_one_EP = 12,
  modred_training = 13,
  modred_saveio_hdf5 = 14,
  Interface = 300,         // Interface timing tasks used for coupling
  Interface_SWEBEM = 301,  // elPaSo/SWEBEM
  Interface_ScaBo = 302,   // elPaSo/ScaBo or Similarf95
  Interface_tBEM = 303,    // elPaSo/tBEM
  initMesh = 706,          // SBFEM timing tasks 706 - 723
  assembleE0E1E2M0 = 707,
  computee1e2m0 = 708,
  initMatrixfield = 709,
  firststep = 710,
  fsInitMatrices = 711,
  fsCompCoeffMat = 712,
  fsSolveRiccati = 713,
  integratingMinf0 = 714,
  schurdecomp = 715,
  nthstep = 716,
  nthstepConvInt = 717,
  nthstep_Lyapunov = 718,
  nthstepTransMinf = 719,
  freeMem = 720,
  dump = 721,
  dump2 = 722
};

//! @brief simple stopwatch to do runtime measurement
//! @author Dirk Clasen
//! @date 26.04.2007
class cStopwatch {
 protected:
  static std::vector<PetscLogDouble>
      m_StartWall;  ///< starting point of each event (wall time)
  static std::vector<PetscLogDouble>
      m_DurationWall;  ///< duration of each event (wall time)

  //! convert duration to a string 'hh:mm:ss'
  void convertToString(char *buffer, eTimingTasks event);

 public:
  //! @brief Constructor
  cStopwatch();
  //! @brief Destructor
  virtual ~cStopwatch();

  //! @brief start to measure time for a specific task
  void startClock(eTimingTasks event);

  //! @brief stop time measurement fo a specific task
  void stopClock(eTimingTasks event);
};

#endif
