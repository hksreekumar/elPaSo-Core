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

#ifndef INFAM_READERH5_H
#define INFAM_READERH5_H

#include <iostream>
#include <vector>

#include "fileh5.h"

#ifdef HAVE_PETSC
#include "petscao.h"
#include "petscconf.h"
#include "petscksp.h"
#include "slepceps.h"
#endif  // HAVE_PETSC

//! @brief H5 class supporting H5 read
//! @author Harikrishnan Sreekumar
//! @date 23.03.2020
//! cOutputH5 reads data from a h5 file.
class cReaderH5 : public cFileH5 {
 public:
  //! @brief Constructor
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  cReaderH5();
  //! @brief Destructor
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  ~cReaderH5();
  //! @brief H5 Routine to read an integer vector
  //! @param _vector Vector handle
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  void readIntegerVector(std::vector<int>& _vector, std::string _nameMainGroup,
                         std::string _nameMatGroup, std::string _nameVecGroup);
  //! @brief H5 Routine to read a dense double matrix
  //! @param _vector Vector handle containing all the matrix details
  //! @param _m Number of rows read
  //! @param _n Number of columns read
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  void readDenseDoubleMatrix(std::vector<double>& _vector, int& _m, int& _n,
                             std::string _nameFirstLevel,
                             std::string _nameSecondLevel,
                             std::string _nameThirdLevel);
  //! @brief H5 Routine to read a dense integer matrix
  //! @param _vector Vector handle containing all the matrix details
  //! @param _m Number of rows read
  //! @param _n Number of columns read
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //! @date 10.06.2021
  //! @author Harikrishnan K. Sreekumar
  void readDenseIntegerMatrix(std::vector<int>& _vector, int& _m, int& _n,
                              std::string _nameFirstLevel,
                              std::string _nameSecondLevel,
                              std::string _nameThirdLevel);

  //! @brief H5 Routine to read an complex matrix
  //! @param _vector Vector handle containing all the matrix details
  //! @param _m Number of rows read
  //! @param _n Number of columns read
  //! @param _nameFirstLevel Main roup name starting with '\'
  //! @param _nameSecondLevel Matrix name
  //! @param _nameThirdLevel Matrix name
  //! @date 09.05.2021
  //! @author Harikrishnan K. Sreekumar
#ifdef PETSC_USE_COMPLEX
  void readDenseComplexMatrix(std::vector<PetscComplex>& _vector, int& _m,
                              int& _n, std::string _nameFirstLevel,
                              std::string _nameSecondLevel,
                              std::string _nameThirdLevel);
#endif

  //! @brief H5 Routine to read an double vector
  //! @param _vector Vector handle
  //! @param _nameMainGroup Main roup name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  void readDoubleVector(std::vector<double>& _vector,
                        std::string _nameMainGroup, std::string _nameMatGroup,
                        std::string _nameVecGroup);

  //! @brief H5 Routine to read an complex vector
  //! @param _vector Vector handle
  //! @param _nameMainGroup Main group name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
#ifdef PETSC_USE_COMPLEX
  void readComplexVector(std::vector<PetscComplex>& _vector,
                         std::string _nameMainGroup, std::string _nameMatGroup,
                         std::string _nameVecGroup);
#endif

  //! @brief H5 Routine to read an complex vector
  //! @param _vector Vector handle
  //! @param _nameMainGroup Main group name starting with '\'
  //! @param _nameMatGroup Matrix name
  //! @param _nameVecGroup Vector name
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  void readNodeCompoundData(std::vector<PetscInt>& _ids,
                            std::vector<double>& _coord,
                            std::string _nameMainGroup,
                            std::string _nameMatGroup,
                            std::string _nameVecGroup);

  //! @brief H5 Routine to read a integer out of an attribute contained in a
  //! group
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @return Read int
  //! @date 25.03.2020
  //! @author Harikrishnan K. Sreekumar
  int readIntegerAttributeFromGroup(std::string _name, std::string _path);

  //! @brief H5 Routine to read a double out of an attribute contained in a
  //! group
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @return Read double
  //! @date 23.04.2021
  //! @author Harikrishnan K. Sreekumar
  double readDoubleAttributeFromGroup(std::string _name, std::string _path);

  //! @brief H5 Routine to read a double out of an attribute contained in a
  //! dataset
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @param _dataset Name of dataset
  //! @return Read double
  //! @date 26.04.2021
  //! @author Harikrishnan K. Sreekumar
  double readDoubleAttributeFromDataset(std::string _name, std::string _path,
                                        std::string _dataset);

  //! @brief H5 Routine to read a string out of an attribute contained in a
  //! group
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @return Read string
  //! @date 25.03.2020
  //! @author Harikrishnan K. Sreekumar
  std::string readStringAttributeFromGroup(std::string _name,
                                           std::string _path);

  //! @brief H5 Routine to read a string out of an attribute contained in a
  //! dataset
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @param _dataset Name of dataset
  //! @return Read string
  //! @date 03.07.2020
  //! @author Harikrishnan K. Sreekumar
  std::string readStringAttributeFromDataset(std::string _name,
                                             std::string _path,
                                             std::string _dataset);

  //! @brief H5 Routine to read a integer out of an attribute contained in a
  //! dataset
  //! @param _name Name of attribute
  //! @param _path Group path to find the attribute
  //! @param _dataset Name of dataset
  //! @return Read integer
  //! @date 20.10.2020
  //! @author Harikrishnan K. Sreekumar
  int readIntegerAttributeFromDataset(std::string _name, std::string _path,
                                      std::string _dataset);

  //! @brief H5 Routine to get number of members of type datasets in given group
  //! path
  //! @param _path Group path
  //! @return number of members
  //! @date 20.10.2020
  //! @author Harikrishnan K. Sreekumar
  int getNumberOfMembersInGroup(std::string _path);

  //! @brief H5 Routine to get number of links in group
  //! @param _path Group path
  //! @return vector of link names
  //! @date 24.01.2022
  //! @author Harikrishnan K. Sreekumar
  PetscInt getNumberOfLinksInGroup(std::string _path);

  //! @brief H5 Routine to read an integer vector from an attribute - dataset
  //! @param _vector Vector handle
  //! @param _path Main roup name starting with '\'
  //! @param _dataset Dataset name
  //! @param _attribute Attribute name
  //! @date 11.05.2021
  //! @author Harikrishnan K. Sreekumar
  void readIntegerVectorFromAttributeDataset(std::vector<int>& _vector,
                                             std::string _path,
                                             std::string _dataset,
                                             std::string _attribute);

  //! @brief H5 Routine to read an integer vector from an attribute - group
  //! @param _vector Vector handle
  //! @param _path Main group name starting with '\'
  //! @param _attribute Attribute name
  //! @date 28.01.2022
  //! @author Harikrishnan K. Sreekumar
  void readIntegerVectorFromAttributeGroup(std::vector<int>& _vector,
                                           std::string _groupname,
                                           std::string _attribute);

  //! @brief H5 Routine to return the name of an object by index
  //! @param _path Group path
  //! @param _index Index of the object
  //! @return name of the member
  //! @date 06.05.2021
  //! @author Harikrishnan K. Sreekumar
  std::string getMemberName(std::string _path, int _index);

 private:
#ifdef HAVE_HDF5
  //! H5 DataSet
  hid_t myDataSet;
  //! H5 DataSpace
  hid_t myDataSpace;
#endif
  //! Complex Struct
  struct elpasoComplexDouble {
    double real;
    double imag;
  };
  //! Nodal Struct
  struct elpasoNodal {
    int Ids;
    double xCoords;
    double yCoords;
    double zCoords;
  };
};
#endif  // !INFAM_READERH5_H
