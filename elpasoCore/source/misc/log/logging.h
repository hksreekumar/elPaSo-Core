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

#ifndef INFAM_LOGGING_H
#define INFAM_LOGGING_H

#include "log.h"

//! @brief class for logging program events
//! @author Dirk Clasen
//! @date 12.12.2006
//! a simple class that is used to write program messages to a text file.
//! only rank 0 is able to write to the logfile; all other ranks to
//! nothing (see documentation of PetscFPrintf() for details).
class cLogging : public virtual cLog {
 private:
 protected:
  static FILE *logfile;  ///< file to which the messages will be written
  PetscErrorCode ierr;   ///< used for error checking

 public:
  cLogging();
  virtual ~cLogging();

  //! @brief write a text to the logfile
  static void trace(const char message[]);

  //! @brief this function can be used like printf() and
  //! calls PetscFPrintf() with the same options
  static void message(const char fmt[], ...);

  //! @brief open file used for logging events
  //! @param filename name of the file
  static void openLogFile(const std::string &filename);

  //! @brief close the file
  static void closeLogFile();

  //! @brief write measured times to log file
  void writeDurations(void);
};

#endif
