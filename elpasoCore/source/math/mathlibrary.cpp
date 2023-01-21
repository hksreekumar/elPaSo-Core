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

#include "mathlibrary.h"
#include <iostream>

#include "../misc/hdf5/outputh5.h"
#include "../misc/hdf5/inputh5.h"
#include "../misc/hdf5/readerh5.h"


#ifdef HAVE_PETSC
void MathLibrary::RetrieveRowVectorsFromMat(Mat& _mat, std::vector<PetscInt>& _ia, std::vector<PetscInt>& _ja, PetscInt shift)
{
    PetscErrorCode ierr; PetscBool done; PetscInt numrows; const PetscInt* ia,* ja;
    ierr = MatGetRowIJ(_mat, shift, PETSC_FALSE, PETSC_FALSE, &numrows, &ia, &ja, &done);
    _ia.insert(_ia.end(), ia, ia + numrows+1);
    _ja.insert(_ja.end(), ja, ja + _ia[_ia.size()-1]-shift);
    MatRestoreRowIJ(_mat, shift, PETSC_FALSE, PETSC_FALSE, &numrows, &ia, &ja, &done);
}

void MathLibrary::RetrieveValueVectorFromMatComplex(Mat& _mat, PetscInt _nnz, std::vector<PetscScalar>& _val)
{
    PetscErrorCode ierr; PetscScalar* val;
    ierr = MatSeqAIJGetArray(_mat, &val);
    _val.insert(_val.end(), val, val + _nnz);
}

void MathLibrary::RetrieveVectorFromVecComplex(Vec& _vec, std::vector<PetscScalar>& _val)
{
    PetscErrorCode ierr; PetscScalar* val; PetscInt size;
    ierr = VecGetArray(_vec, &val);
    ierr = VecGetSize(_vec, &size);
    _val.insert(_val.end(), val, val + size);
}

void MathLibrary::ExportPetscMatToHDF5Complex(Mat& _mat, std::string _systemtype, std::string _filename)
{
#ifdef HAVE_HDF5
    // Intitialize
    std::vector<PetscInt> ia, ja;
    std::vector<PetscScalar> val;
    PetscInt m, n;

    // Retrieve Information of Petsc Mat
    RetrieveRowVectorsFromMat(_mat, ia, ja, 1);
    RetrieveValueVectorFromMatComplex(_mat, ja.size(), val);
    MatGetSize(_mat, &m, &n);

    // Export to HDF5
    cOutputH5 systemoutput_hdf5;
    systemoutput_hdf5.setGlobalFileName(_filename);

    if (systemoutput_hdf5.exists_hdf5())
        systemoutput_hdf5.openContainer(ELPASO_H5_READWRITE);
    else
        systemoutput_hdf5.openContainer(ELPASO_H5_READWRITE_FORCE);

    systemoutput_hdf5.createGroup("/SystemMatrices");
    systemoutput_hdf5.createSecoGroup("/SystemMatrices", _systemtype);
    systemoutput_hdf5.writeIntegerAttributeToGroup("numrow", "/SystemMatrices" + _systemtype, m);
    systemoutput_hdf5.writeIntegerAttributeToGroup("numcol", "/SystemMatrices" + _systemtype, n);
    systemoutput_hdf5.appendIntegerVector(ia, "/SystemMatrices", _systemtype, "/vecCsrIa");
    systemoutput_hdf5.appendIntegerVector(ja, "/SystemMatrices", _systemtype, "/vecCsrJa");
#ifdef PETSC_USE_COMPLEX
    systemoutput_hdf5.appendComplexVector(val, "/SystemMatrices", _systemtype, "/cmpCsrVal");
#endif
    systemoutput_hdf5.closeContainer();
#endif // HAVE_HDF5
}


/*BEGIN_NO_COVERAGE*/
void MathLibrary::ExportPetscVecToHDF5Complex(Vec& _vec, std::string _namePrimGroup, std::string _nameSecoGroup, std::string _nameTerGroup, std::string _filename)
{
    // Initialize
    std::vector<PetscScalar> vec;

    // Retrieve Information of Petsc Vec
    RetrieveVectorFromVecComplex(_vec, vec);
    // Export to HDF5
    cOutputH5 systemoutput_hdf5;
    systemoutput_hdf5.setGlobalFileName(_filename);

    if (systemoutput_hdf5.exists_hdf5())
        systemoutput_hdf5.openContainer(ELPASO_H5_READWRITE);
    else
        systemoutput_hdf5.openContainer(ELPASO_H5_READWRITE_FORCE);
    
    systemoutput_hdf5.createGroup(_namePrimGroup);
    systemoutput_hdf5.createSecoGroup(_namePrimGroup, _nameSecoGroup);
#ifdef PETSC_USE_COMPLEX
    systemoutput_hdf5.appendComplexVector(vec, _namePrimGroup, _nameSecoGroup, _nameTerGroup);
#endif
    systemoutput_hdf5.closeContainer();
}
/*END_NO_COVERAGE*/

#endif