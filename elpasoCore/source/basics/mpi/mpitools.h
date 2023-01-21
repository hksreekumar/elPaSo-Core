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

#ifndef CMPITOOLS_H
#define CMPITOOLS_H

#include <mpi.h>

#include <sstream>
#include <vector>

#include "../exceptions/exceptions.h"

//! @brief functions to make mpi communications a lot easier
//! @author Marco Schauer
//! @note thx to SirAnn ;)
//! @date 28.04.2009
namespace cMpiTools {
//! @brief send a single value "value" of MPI_Datatype
template <class T>
inline void sendSingleValue(T& value, MPI_Datatype datatype, int dest, int tag,
                            MPI_Comm comm);

//! @brief receives a single value "value" of MPI_Datatype
template <class T>
inline void receiveSingleValue(T& value, MPI_Datatype datatype, int source,
                               int tag, MPI_Comm comm);

//! @brief receives and returns a single value of MPI_Datatype
//! expample: int value =
//! PbMpiTools::receiveSingleValue<int>(MPI::INT,0,10,comm);
template <class T>
inline T receiveSingleValue(MPI_Datatype datatype, int source, int tag,
                            MPI_Comm comm);

//! @brief sends bool value (doesn't work with template, why ever... stupid MPI)
inline void sendBoolValue(bool value, int dest, int tag, MPI_Comm comm);

//! @brief receives bool value (doesn't work with template, why ever... stupid
//! MPI)
inline bool receiveBoolValue(int source, int tag, MPI_Comm comm);

//! @brief sends bool value (doesn't work with template, why ever... stupid MPI)
inline void sendStringValue(const std::string& value, int dest, int tag,
                            MPI_Comm comm);

//! @brief receives bool value (doesn't work with template, why ever... stupid
//! MPI)
inline std::string receiveStringValue(int source, int tag, MPI_Comm comm);

//! @brief send a vector of MPI_Datatype
template <class T>
inline void sendVector(std::vector<T>& v, MPI_Datatype datatype, int dest,
                       int tag, MPI_Comm comm);

//! @brief receive a std::vector of MPI_Datatype
template <class T>
inline void receiveVector(std::vector<T>& v, MPI_Datatype datatype, int source,
                          int tag, MPI_Comm comm);

//! @brief receive a vector of MPI_Datatype and adds this vector to existing
//! vector ans returns number of received elements
template <class T>
inline int receiveVectorAndAddToVector(std::vector<T>& v, MPI_Datatype datatype,
                                       int source, int tag, MPI_Comm comm);

//! @brief send a std::vector of strings
inline void sendStringVector(const std::vector<std::string>& v, int dest,
                             int tag, MPI_Comm comm);

//! @brief send a vector of strings
inline void receiveStringVector(std::vector<std::string>& v, int dest, int tag,
                                MPI_Comm comm);
};  // namespace cMpiTools

//! @brief send a single value of MPI_Datatype
template <class T>
void cMpiTools::sendSingleValue(T& value, MPI_Datatype datatype, int dest,
                                int tag, MPI_Comm comm) {
  MPI_Send(&value, 1, datatype, dest, tag, comm);
}

template <class T>
void cMpiTools::receiveSingleValue(T& value, MPI_Datatype datatype, int source,
                                   int tag, MPI_Comm comm) {
  MPI_Status status;
  MPI_Recv(&value, 1, datatype, source, tag, comm, &status);
}

template <class T>
T cMpiTools::receiveSingleValue(MPI_Datatype datatype, int source, int tag,
                                MPI_Comm comm) {
  MPI_Status status;
  T value;
  MPI_Recv(&value, 1, datatype, source, tag, comm, &status);
  return value;
}

//! @brief send a bool value (bool doesn't work with template, why ever)
void cMpiTools::sendBoolValue(bool value, int dest, int tag, MPI_Comm comm) {
  short dummy;
  if (value)
    dummy = 1;
  else
    dummy = 0;

  MPI_Send(&dummy, 1, MPI_SHORT, dest, tag, comm);
}

bool cMpiTools::receiveBoolValue(int source, int tag, MPI_Comm comm) {
  MPI_Status status;
  short dummy;
  MPI_Recv(&dummy, 1, MPI_SHORT, source, tag, comm, &status);
  return (dummy == 1);
}

//! @brief sends bool value (doesn't work with template, why ever... stupid MPI)
void cMpiTools::sendStringValue(const std::string& value, int dest, int tag,
                                MPI_Comm comm) {
  std::vector<char> vec;
  for (std::size_t i = 0; i < value.size(); i++) vec.push_back(value[i]);
  cMpiTools::sendVector(vec, MPI_CHAR, dest, tag, comm);
}

//! @brief receives bool value (doesn't work with template, why ever... stupid
//! MPI)
std::string cMpiTools::receiveStringValue(int source, int tag, MPI_Comm comm) {
  std::vector<char> vec;
  cMpiTools::receiveVector(vec, MPI_CHAR, source, tag, comm);
  std::string str;
  for (std::size_t i = 0; i < vec.size(); i++) str += vec[i];

  return str;
}

//! @brief send a vector of MPI_Datatype
template <class T>
void cMpiTools::sendVector(std::vector<T>& v, MPI_Datatype datatype, int dest,
                           int tag, MPI_Comm comm) {
  // send size
  int size = (int)v.size();

  MPI_Send(&size, 1, MPI_INT, dest, tag, comm);

  if (size > 0) MPI_Send(&v[0], size, datatype, dest, tag, comm);
}

//! @brief receive a vector of MPI_Datatype
template <class T>
void cMpiTools::receiveVector(std::vector<T>& v, MPI_Datatype datatype,
                              int source, int tag, MPI_Comm comm) {
  MPI_Status status;
  int size = 0;
  MPI_Recv(&size, 1, MPI_INT, source, tag, comm, &status);
  v.resize(size);
  if (size > 0) MPI_Recv(&v[0], size, datatype, source, tag, comm, &status);
}

//! @brief receive a vector of MPI_Datatype and adds this vector to existing
//! vector
//! @return value is size of received elements
template <class T>
int cMpiTools::receiveVectorAndAddToVector(std::vector<T>& v,
                                           MPI_Datatype datatype, int source,
                                           int tag, MPI_Comm comm) {
  MPI_Status status;
  int incommingSize;
  MPI_Recv(&incommingSize, 1, MPI_INT, source, tag, comm, &status);
  int oldSize = (int)v.size();
  v.resize(oldSize + incommingSize);
  if (incommingSize > 0)
    MPI_Recv(&v[oldSize], incommingSize, datatype, source, tag, comm, &status);
  return incommingSize;
}

//! @brief send a vector of strings
void cMpiTools::sendStringVector(const std::vector<std::string>& v, int dest,
                                 int tag, MPI_Comm comm) {
  // send size
  int stringVectorSize = (int)v.size();
  MPI_Send(&stringVectorSize, 1, MPI_INT, dest, tag, comm);
  if (stringVectorSize > 0) {
    std::vector<int> singleStringSizes(stringVectorSize + 1);
    int nofChars = 0;
    for (int i = 0; i < stringVectorSize; i++)
      nofChars += singleStringSizes[i] = (int)v[i].length();

    singleStringSizes[stringVectorSize] = nofChars;
    MPI_Send(&singleStringSizes[0], stringVectorSize + 1, MPI_INT, dest, tag,
             comm);

    std::vector<char> charVector(nofChars);
    int pos = 0;
    for (int i = 0; i < stringVectorSize; i++)
      for (int j = 0; j < singleStringSizes[i]; j++)
        charVector[pos++] = v[i][j];

    MPI_Send(&charVector[0], nofChars, MPI_CHAR, dest, tag, comm);
  }
}

//! @brief send a vector of strings
void cMpiTools::receiveStringVector(std::vector<std::string>& v, int source,
                                    int tag, MPI_Comm comm) {
  MPI_Status status;
  // send size
  int stringVectorSize;
  MPI_Recv(&stringVectorSize, 1, MPI_INT, source, tag, comm, &status);

  v.clear();
  v.resize(stringVectorSize);

  if (stringVectorSize > 0) {
    std::vector<int> singleStringSizes(stringVectorSize + 1);

    MPI_Recv(&singleStringSizes[0], stringVectorSize + 1, MPI_INT, source, tag,
             comm, &status);

    int nofChars = singleStringSizes[stringVectorSize];
    std::vector<char> charVector(nofChars);

    MPI_Recv(&charVector[0], nofChars, MPI_CHAR, source, tag, comm, &status);

    int pos = 0;
    for (int i = 0; i < stringVectorSize; i++)
      for (int j = 0; j < singleStringSizes[i]; j++)
        v[i].push_back(charVector[pos++]);
  }
}

#endif  // CMPITOOLS_H
