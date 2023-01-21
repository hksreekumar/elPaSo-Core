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

#ifndef INFAM_OUTPUTH5_H
#define INFAM_OUTPUTH5_H

//#include "outputfile.h"
#include <vector>

#include "fileh5.h"

#ifdef HAVE_PETSC
#include "petscao.h"
#include "petscconf.h"
#include "petscksp.h"
#include "slepceps.h"

#endif  // HAVE_PETSC

//! @brief H5 class supporting H5 write
//! @author Harikrishnan Sreekumar
//! @date 18.03.2020
//! cOutputH5 writes data to a h5 file.

class cOutputH5 : public cFileH5 {
 public:
  //! @brief H5 Routine to create a main group
  //! @param _nameMainGroup Main group name starting with '\'
  //! @date 19.03.2020
  //! @author Harikrishnan K. Sreekumar
  void createGroup(std::string _nameMainGroup);
  //! @brief H5 Routine to create a secondary group
  //! @param _nameMainGroup Main group name starting with '\'
  //! @param _nameSecoGroup Secondary group name starting with '\'
  //! @date 19.03.2020
  //! @author Harikrishnan K. Sreekumar
  void createSecoGroup(std::string _nameMainGroup, std::string _nameSecoGroup);
  //! @brief H5 Routine to append an integer vector
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //!      The structure of resulting H5 :
  //!           H5 RootFile:
  //!                   -- _nameMainGroup
  //!                          -- _nameMatGroup
  //!                                   -- _nameVecGroup
  //! @date 19.03.2020
  //! @author Harikrishnan K. Sreekumar
  void appendIntegerVector(std::vector<int>& _vector,
                           std::string _nameMainGroup,
                           std::string _nameMatGroup,
                           std::string _nameVecGroup);

  //! @brief H5 Routine to write an integer vector as an attribute
  //! @param _vector Vector handle
  //! @param _path Path to group
  //! @param _attribute Attribute name
  //! @date 28.01.2022
  //! @author Harikrishnan K. Sreekumar
  //! @see testOutputHDF5.h
  void writeIntegerVectorAttributeToGroup(std::vector<int>& _vector,
                                          std::string _pathGroup,
                                          std::string _attribute);

  //! @brief H5 Routine to append an complex vector
  //! @param _vector Vector to be appendend to HDF5 file
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //!      The structure of resulting H5 :
  //!           H5 RootFile:
  //!                   -- _nameMainGroup
  //!                          -- _nameMatGroup
  //!                                   -- _nameVecGroup
  //! @date 19.03.2020
  //! @author Harikrishnan K. Sreekumar
#ifdef PETSC_USE_COMPLEX
  void appendComplexVector(std::vector<PetscComplex>& _vector,
                           std::string _nameMainGroup,
                           std::string _nameMatGroup,
                           std::string _nameVecGroup);
#endif
  //! @brief H5 Routine to append an complex vector using hyperslabs (parallel)
  //! @param _vector Vector to be appendend to HDF5 file
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //!      The structure of resulting H5 :
  //!           H5 RootFile:
  //!                   -- _nameMainGroup
  //!                          -- _nameMatGroup
  //!                                   -- _nameVecGroup
  //! @param _start Start index
  //! @param _end   Ending index
  //! @date 22.09.2021
  //! @author Harikrishnan K. Sreekumar
#ifdef PETSC_USE_COMPLEX
  void appendComplexVector(std::vector<PetscComplex>& _vector,
                           std::string _nameMainGroup,
                           std::string _nameMatGroup, std::string _nameVecGroup,
                           PetscInt _start, PetscInt _end);
#endif
  //! @brief H5 Routine to append a complex matrix
  //! @param _vector Vector of matrix
  //! @param _m Number of rows
  //! @param _n Number of columns
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Matrix name
  //!      The structure of resulting H5 :
  //!           H5 RootFile:
  //!                   -- _nameMainGroup
  //!                          -- _nameMatGroup
  //!                                   -- _nameVecGroup
  //! @param _isRowMajor True if row major, False if column major
  //! @date 25.08.2020
  //! @author Harikrishnan K. Sreekumar
  void appendComplexMatrix(std::vector<PetscScalar>& _vectorMatrix, size_t _m,
                           size_t _n, std::string _nameMainGroup,
                           std::string _nameMatGroup, std::string _nameVecGroup,
                           PetscBool _isRowMajor);
  //! @brief H5 Routine to append a dense real matrx
  //! @param _mat Matrix to be exported
  //! @param _m Number of rows
  //! @param _n Number of columns
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //!      The structure of resulting H5 :
  //!           H5 RootFile:
  //!                   -- _nameMainGroup
  //!                          -- _nameMatGroup
  //!                                   -- _nameVecGroup
  //! @date 25.04.2020
  //! @author Harikrishnan K. Sreekumar
#ifdef HAVE_PETSC
  void appendRealDenseMatrix(std::vector<PetscReal>& _matrix, int _m, int _n,
                             std::string _nameMainGroup,
                             std::string _nameMatGroup,
                             std::string _nameVecGroup);
#endif
  //! @brief H5 Routine to append a dense integer matrx
  //! @param _mat Matrix to be exported
  //! @param _m Number of rows
  //! @param _n Number of columns
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //!      The structure of resulting H5 :
  //!           H5 RootFile:
  //!                   -- _nameMainGroup
  //!                          -- _nameMatGroup
  //!                                   -- _nameVecGroup
  //! @date 14.12.2021
  //! @author Harikrishnan K. Sreekumar
#ifdef HAVE_PETSC
  void appendIntegerDenseMatrix(std::vector<PetscInt>& _matrix, int _m, int _n,
                                std::string _nameMainGroup,
                                std::string _nameMatGroup,
                                std::string _nameVecGroup);
#endif
  //! @brief Constructor
  //! @date 19.03.2020
  //! @author Harikrishnan K. Sreekumar
  cOutputH5();
  //! @brief Destructor
  //! @date 19.03.2020
  //! @author Harikrishnan K. Sreekumar
  ~cOutputH5();
  //! @brief Function to return the singleton instance
  //! @date 24.03.2020
  //! @author Harikrishnan K. Sreekumar
  /*static cOutputH5* getInstance() {
          if (!mySingleton) {
              mySingleton = new cOutputH5();
          }
          return mySingleton;
  }*/

  //! @brief H5 Routine to write a string to an attribute
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @param _value Value of attribute
  //! @date 25.03.2020
  //! @author Harikrishnan K. Sreekumar
  void writeStringAttributeToGroup(std::string _name, std::string _path,
                                   std::string _value);
  //! @brief H5 Routine to write a integer to an attribute
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @param _value Value of attribute
  //! @date 25.05.2021
  //! @author Harikrishnan K. Sreekumar
  void writeIntegerAttributeToGroup(std::string _name, std::string _path,
                                    int _value);
  //! @brief H5 Routine to write a integer to a dataset
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @param _value Value of attribute
  //! @param _dataset Name of dataset
  //! @date 21.01.2022
  //! @author Harikrishnan K. Sreekumar
  void writeIntegerAttributeToDataset(std::string _name, std::string _path,
                                      int _value, std::string _dataset);

 private:
#ifdef HAVE_HDF5
  //! H5 Group
  hid_t myGroup;
  //! H5 DataSet
  hid_t myDataSet;
  //! H5 DataSpace
  hid_t myDataSpace;
#endif
  struct elpasoComplexDouble {
    double real;
    double imag;
  };
};
#endif