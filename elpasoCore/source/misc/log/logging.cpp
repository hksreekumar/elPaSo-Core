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

#include "logging.h"

#ifndef NO_BUILDINFO_H
#include "buildinfo.h"
#endif

cLogging::cLogging()
{
  // empty
}

cLogging::~cLogging()
{
  // empty
}

// --- initialise static members
FILE* cLogging::logfile = NULL;

void cLogging::trace(const char message[])
{
#ifndef ELPASO_TEST
  PetscFPrintf(PETSC_COMM_WORLD, logfile, "%s\n", message);
#endif // ELPASO_TEST
}


void cLogging::message(const char fmt[], ...)
{
#ifndef ELPASO_TEST
  char buffer[1024];

  va_list args;
  va_start(args, fmt);
  vsprintf(buffer, fmt, args);
  va_end(args);

  PetscFPrintf(PETSC_COMM_WORLD, logfile, buffer);
#endif // ELPASO_TEST
}

void cLogging::openLogFile(const std::string &filename)
{
  if (logfile != NULL) {
    trace("*** YOU TRIED TO REOPEN THE LOGFILE!! CHECK YOUR CODE");
    return;
  }

  PetscFOpen(PETSC_COMM_WORLD, filename.c_str(), "w", &logfile);

  trace  ( "========================================================================" );
  trace  ( "  ELPASO " );
  trace  ( "------------------------------------------------------------------------" );

  message( "  Version         : %s\n", GIT_COUNT);
  message( "  Git-Hash        : %s\n", GIT_COMMIT_HASH );
  message( "  Compiled  Date  : %s\n", __DATE__ );
  message( "            Time  : %s\n", __TIME__ );
  trace  ( "------------------------------------------------------------------------" );
  message( "  Versions of core packages: \n");
#ifdef OMPI_MAJOR_VERSION
  message( "  OpenMPI_VER     : %d.%d.%d \n", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION, OMPI_RELEASE_VERSION );
#endif
#ifdef I_MPI_VERSION
  message( "  IntelMPI_VER    : %s \n", I_MPI_VERSION );
#endif
  message( "  PETSC_ARCH      : %s\n", PETSC_ARCH ); //Since PETSC 3.1 we have to use PETSC_ARCH. Previous PETSC versions need PETSC_ARCH_NAME
  message( "  PETSC_VER       : %d.%d.%d-p%d \n", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR, PETSC_VERSION_PATCH );
  message( "  SLEPC_VER       : %d.%d.%d-p%d \n", SLEPC_VERSION_MAJOR, SLEPC_VERSION_MINOR, SLEPC_VERSION_SUBMINOR, SLEPC_VERSION_PATCH );
#ifdef HDF5_VERSION
  message( "  HDF5_VER        : %s \n", HDF5_VERSION );
#endif
  message( "  ARPACK/PARPACK  : 2.1 (Release: 1995-11-15) \n" ); /* Can be found in version.h of ARPACK/SRC/ directory */
  message( "  *** CAUTION: ARPACK/PARPACK  version cannot be properly tracked! *** \n" );
  trace  ( "========================================================================" );
}

void cLogging::closeLogFile()
{
  PetscFClose(PETSC_COMM_WORLD, logfile);
}

void cLogging::writeDurations(void)
{
  char buffer[128];

  trace  ( "------------------------------------------------------------------------" );
  message( "  RESOURCE USAGE:         [dd:hh:mm:ss] \n");
  convertToString(buffer, Overall);          message("    Overall ............: %s\n", buffer);
  convertToString(buffer, reorder);          message("     reorder ...........: %s\n", buffer);
  convertToString(buffer, SetupKM);          message("     compute K and M ...: %s\n", buffer);
  convertToString(buffer, Solver);           message("     solve system ......: %s\n", buffer);
  ierr =  PetscPrintf(PETSC_COMM_WORLD,"Solution time: %s\n", buffer);
  convertToString(buffer, Interface);        message("     interface .........: %s\n", buffer);
  convertToString(buffer, Interface_SWEBEM); message("      interface_SWEBEM .: %s\n", buffer);
  convertToString(buffer, Interface_ScaBo);  message("      interface_ScaBo ..: %s\n", buffer);
  convertToString(buffer, initMatrixfield);  message("        SBFEM M^infty ..: %s\n", buffer);
  convertToString(buffer, Interface_tBEM);   message("      interface_tBEM ...: %s\n", buffer);
  convertToString(buffer, output);           message("     output ............: %s\n", buffer);
  trace("-- MODRED --------------------------------------------------------------");
  convertToString(buffer, modred_overall);   message("    Overall ............: %s\n", buffer);
  convertToString(buffer, modred_training);  message("     training ..........: %s\n", buffer);
  convertToString(buffer, modred_saveio_hdf5);   message("     io mem save .......: %s\n", buffer);
  //ierr = PetscPrintf(PETSC_COMM_WORLD, "Modred time: %s\n", buffer);
  trace  ( "------------------------------------------------------------------------" );
}
