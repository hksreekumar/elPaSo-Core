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

#ifndef INFAM_FILEH5_H
#define INFAM_FILEH5_H

#include <sys/stat.h>

#include <string>

#ifdef HAVE_HDF5
#include "hdf5.h"
#endif  // HAVE_HDF5

//! @brief HDF5 file read mode
enum ELPASO_H5MODE {
  ELPASO_H5_READONLY,
  ELPASO_H5_READWRITE,
  ELPASO_H5_READWRITE_FORCE
};

//! @brief Implements pattern for the global elpaso H5 instance
//! @author Harikrishnan Sreekumar
//! @date 23.03.2020
//! cFileH5 is the base class for input and output routines of H5
class cFileH5 {
 public:
  //! @brief Constructor
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  cFileH5();
  //! @brief Destructor
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  ~cFileH5();
  //! @brief Base H5 Routine to set h5 file name
  //! @param _filename Name of file
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  void setGlobalFileName(std::string _filename) { myFileName = _filename; }
  //! @brief Base H5 Routine to get h5 file name
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  std::string getGlobalFileName() { return myFileName; }
  //! @brief Base H5 Routine to close a h5 file
  //! @param _mode Mode for opening container [READ-ONLY or READ/WRITE]
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  void openContainer(ELPASO_H5MODE _mode);
  //! @brief Base H5 Routine to close a h5 file
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  virtual void closeContainer();

  //! @brief Function to return the static filename instance
  //! @date 24.03.2020
  //! @author Harikrishnan K. Sreekumar
#ifdef HAVE_HDF5
  hid_t getFileInstance() {
    if (my_file_id) {
      return my_file_id;
    }
    return NULL;
  };
#endif  // HAVE_HDF5

  inline bool exists_hdf5() {
    const std::string name = getGlobalFileName();
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
  }

 private:
#ifdef HAVE_HDF5
  // file access propery
  hid_t my_plist_id;  // property list identifier

  /// H5 Root file
  hid_t my_file_id;

#endif
  /// File name
  std::string myFileName;
};
#endif