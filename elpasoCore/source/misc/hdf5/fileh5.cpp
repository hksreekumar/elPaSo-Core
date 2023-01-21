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

#include "fileh5.h"

#include <iostream>


cFileH5::cFileH5() {
  // Empty
}

cFileH5::~cFileH5() {}

void cFileH5::openContainer(ELPASO_H5MODE _mode) {
#ifdef HAVE_HDF5
  my_plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(my_plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  hid_t access_plist;

  switch (_mode) {
    case ELPASO_H5_READONLY:
      // Existing file is opened with read-only access. If file
      // does not exist, H5Fopen fails.
      my_file_id = H5Fopen(myFileName.c_str(), H5F_ACC_RDONLY, my_plist_id);
      break;
    case ELPASO_H5_READWRITE:
      // Existing file is opened with read - write access. If file
      // does not exist, H5Fopen fails.
      my_file_id = H5Fopen(myFileName.c_str(), H5F_ACC_RDWR, my_plist_id);
      break;
    case ELPASO_H5_READWRITE_FORCE:
      // If file already exists, file is opened with read-write
      // access and new data overwrites existing data,
      // destroying all prior content, i.e., file content is truncated
      // upon opening. If file does not exist, it is created and
      // opened with read-write access.
      my_file_id = H5Fcreate(myFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                             my_plist_id);
      break;
    default:
      std::cerr << "	! Invalid H5 open mode !" << std::endl;
      break;
  }
  H5Pclose(my_plist_id);
#endif  // HAVE_HDF5
}

void cFileH5::closeContainer() {
#ifdef HAVE_HDF5
  H5Fclose(my_file_id);
#endif
}
