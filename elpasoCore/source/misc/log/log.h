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

#ifndef INFAM_LOG_H
#define INFAM_LOG_H

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "../../basics/watch/stopwatch.h"
#include "../counter.h"
#include "mpi.h" /* OpenMPI version variables for logging.cpp */
#include "petsc.h"
#include "petscversion.h" /* Petsc version variables for logging.cpp */
#include "slepcversion.h" /* Slepc version variables for logging.cpp */

//! @brief class for logging some data or messages
//! @author Marco Schauer
//! @date 01.10.2008
//! a simple class that is used to write some data or messages to a text file.
//! only rank 0 is able to write to the logfile; all other ranks to
//! nothing (see documentation of PetscFPrintf() for details).
class cLog : public virtual cStopwatch {
 private:
  bool m_mute;

 protected:
 public:
  cLog();
  virtual ~cLog();
  void setMute(bool mute) { m_mute = mute; };
  bool getMute() { return m_mute; };
};

#endif
