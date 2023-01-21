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

#include "femparserfactory.h"
#include "femparserinterface.h"
#include "femparserhdf5.h"

cFemParserFactory::cFemParserFactory()
{
    // empty
}

cFemParserFactory::~cFemParserFactory()
{
    // empty
}

cFemParserInterface* cFemParserFactory::createParser(std::string _fileExtension)
{
    // create instance
    cFemParserInterface* ptr = nullptr;
    if (_fileExtension == "hdf5"){
        PetscPrintf(PETSC_COMM_SELF, "Using hdf5 parser\n");
        ptr = new cFemParserHDF5;
    }
    else
        PetscPrintf(PETSC_COMM_SELF, "Unidentified parser\n");
        
    return ptr;
}
