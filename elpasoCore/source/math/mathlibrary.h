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

#ifdef HAVE_PETSC
#include "petscao.h"
#include "petscconf.h"
#include "petscksp.h"
#include "slepceps.h"

#endif
#include <vector>

#ifdef HAVE_INTELMKL
#include "mkl.h"
#endif

//! @brief function definitions for all useful math functions
//! @date 28.10.2019
//! @author Harikrishnan K. Sreekumar
namespace MathLibrary {
#ifdef HAVE_PETSC
    /**
    * @brief Petsc routine to retrieve CSR ia and ja vectors from mat
    * @param _mat Matrix to be initialized
    * @param _ia IA vector of the CSR format
    * @param _ja JA vector of the CSR format
    * @param _n Number of rows
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    * @see testmathlibrary_matsparse.h, testmathlibrary_matdense.h
    **/
  void RetrieveRowVectorsFromMat(Mat& _mat, std::vector<PetscInt>& _ia, std::vector<PetscInt>& _ja, PetscInt shift);
  
    /**
    * @brief Petsc routine to retrieve CSR value vector from mat
    * @param _mat Matrix to be initialized
    * @param _nnz Number of non-zeroes
    * @param _val Value vector of the CSR format [Complex Valued]
    * @date 19.03.2020
    * @author Harikrishnan K. Sreekumar
    * @see testmathlibrary_matsparse.h, testmathlibrary_matdense.h
    **/
    void RetrieveValueVectorFromMatComplex(Mat& _mat, PetscInt _nnz, std::vector<PetscScalar>& _val);

    /**
    * @brief Petsc routine to retrieve vector from PETSc Vec
    * @param _vec Vector
    * @param _val Vector retrieved [Complex Valued]
    * @date 06.04.2020
    * @author Harikrishnan K. Sreekumar
    * @see testmathlibrary_matsparse.h, testmathlibrary_matdense.h
    **/
    void RetrieveVectorFromVecComplex(Vec& _vec, std::vector<PetscScalar>& _val);

    /**
    * @brief Exports Petsc MAT to an H5 File in CSR format
    * @param _mat Matrix to be exported
    * @param _systemtype H5 System type
    * @param _filename HDF5 file name
    * @date 23.03.2020
    * @author Harikrishnan K. Sreekumar
    * @see testmathlibrary_matsparse.h, testmathlibrary_matdense.h
    **/
    void ExportPetscMatToHDF5Complex(Mat& _mat, std::string _systemtype, std::string _filename);

    /**
    * @brief Exports Petsc VEC to an H5 File 
    * @param _vec Vector to be exported
    * @param _namePrimGroup Primary group
    * @param _nameSecoGroup Secondary group
    * @param _nameTerGroup Tertiary group
    * @param _filename HDF5 file name
    * @date 06.04.2020
    * @author Harikrishnan K. Sreekumar
    **/
    void ExportPetscVecToHDF5Complex(Vec& _vec, std::string _namePrimGroup, std::string _nameSecoGroup, std::string _nameTerGroup, std::string _filename);
#endif
}

//! this is a wrapper function of PetscError. PETSc specifies CHKERRQ but this
//! macro contains a return value. Therefore it is useless within functions
//! specified as void. To get rid of the problem this wrapper is used throughout
//! the whole code. So never use CHKERRQ - you'll get funny compiler messages.
#ifdef PETSC_USE_ERRORCHECKING
inline void INFAMCHKERRQ(PetscErrorCode a) {
  if (a) {
    PetscError(MPI_COMM_WORLD, __LINE__, __FUNCT__, __FILE__, a,
               PETSC_ERROR_INITIAL, "**** OOPS **** ");
    MPI_Abort(MPI_COMM_WORLD, a);
  };
}
#else
#define INFAMCHKERRQ(a) {};
#endif
