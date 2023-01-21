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

#ifndef INFAM_INPUTH5_H
#define INFAM_INPUTH5_H

#include <vector>

#include "readerh5.h"

#ifdef HAVE_PETSC
#include "petscao.h"
#include "petscconf.h"
#include "petscksp.h"
#include "slepceps.h"
#endif  // HAVE_PETSC

//! @brief Singleton H5 class supporting H5 read
//! @author Harikrishnan Sreekumar
//! @date 23.03.2020
//! cOutputH5 reads data from a h5 file.

class cInputSingletonH5 : public cReaderH5 {
 public:
  //! @brief Constructor
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  cInputSingletonH5();
  //! @brief Destructor
  //! @date 23.03.2020
  //! @author Harikrishnan K. Sreekumar
  ~cInputSingletonH5();
  //! @brief Function to return the singleton instance
  //! @date 19.03.2020
  //! @author Harikrishnan K. Sreekumar
  static cInputSingletonH5* getInstance() {
    if (!mySingleton) {
      mySingleton = new cInputSingletonH5();
    }
    return mySingleton;
  }

 private:
  //! Singleton
  static cInputSingletonH5* mySingleton;
};
#endif  // !INFAM_INPUTH5_H
