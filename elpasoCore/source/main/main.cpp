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

#include "../misc/parser/femparserfactory.h"
#include "../misc/parser/femparserinterface.h"

#include <iostream>

const char help[] = "\n\
  **********************************************************\n\
  elPaSo: elementary Parallel Solver\n\
  Vibroacoustic Finite Element Code\n\
  Technische Universitaet Braunschweig\n\
  Institute for Acoustics\n\
  Langer Kamp 19\n\
  D-38106 Braunschweig, Germany\n\
  **********************************************************\n\n";

/**
* Main program of the finite element code
*/
int main(int argc, char* argv[])
{
  // --- initialize MPI, PETSc, and SLEPc
  SlepcInitialize(&argc,&argv,0,help);
  PetscMemorySetGetMaximumUsage();

  char                  cFilename[100];
  std::string           myFilename, myFileExtension;
  std::string           logname;       // filename of logfile
  PetscBool             found = PETSC_FALSE;
  cFemParserFactory     myParserFactory;
  cFemParserInterface*  myParser = 0;
  cProblem              myProblem;     // stores data related to current problem (mesh, analysis parameters)
  cLogging              Logging;
  PetscBool             optCompute = PETSC_FALSE;

  // -----------------------------------------------------------------------
  //   check, if an inputfile is given. Otherwise stop program
  // -----------------------------------------------------------------------
  PetscOptionsGetString(PETSC_NULL,"", "-inp", cFilename, 100, &found);
  if (found == PETSC_FALSE) {
    PetscPrintf(PETSC_COMM_WORLD, "\n\nNO INPUTFILE - please specify inputfile with -inp option. GOOD BYE ...\n\n");
    SlepcFinalize();
    return 1;
  }

  // -----------------------------------------------------------------------
  //   strip suffix from name of inputfile
  // -----------------------------------------------------------------------
  myFilename = std::string( cFilename );
  myFileExtension = myFilename.substr(myFilename.find_last_of(".") + 1) ;
  if (myFilename.rfind(myFileExtension) == (myFilename.length()-myFileExtension.length()))
    myFilename.erase(myFilename.length()-myFileExtension.length()-1);

  // -----------------------------------------------------------------------
  //   establish parser
  // -----------------------------------------------------------------------
  myParser = myParserFactory.createParser(myFileExtension);
  if (myParser == nullptr)
  {
    PetscPrintf(PETSC_COMM_WORLD, "\n\nelPaSo cannot identify the current file-format. GOOD BYE ...\n\n");
    SlepcFinalize();
    return 1;
  }

  // -----------------------------------------------------------------------
  //   open logfile
  // -----------------------------------------------------------------------
  Logging.openLogFile( myFilename + ".log.0" );

  // -----------------------------------------------------------------------
  //   actual parsing
  // -----------------------------------------------------------------------
  myProblem.setFilename(myFilename + "." + myFileExtension);
  myParser->parseInputFile(myProblem);
  myProblem.getMesh()->checkElementSizes();
  myProblem.getAnalysis()->setBaseName(myFilename);

  // ------------------------------------------------------------------------
  //   perform computation on the given mesh
  // ------------------------------------------------------------------------
  PetscOptionsHasName(PETSC_NULL,PETSC_NULL, "-c", &optCompute);
  if (optCompute == PETSC_TRUE) //user added option "-c"
  {
      // --- perform the requested computations
      myProblem.getAnalysis()->FullRun(*(myProblem.getMesh()));
  }
  
  Logging.stopClock(Overall);
  Logging.writeDurations();

  // -----------------------------------------------------------------------
  //   deallocations
  // -----------------------------------------------------------------------
  myParser->closeInputFile();
  myParser = NULL;
  delete myParser;

  // -- close log files, shutdown MPI
  Logging.closeLogFile();
  SlepcFinalize();
  return 0;
}