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

#ifndef INFAM_EXCEPTIONS_H
#define INFAM_EXCEPTIONS_H

#include <sstream>
#include <string>

//! @brief handling exceptions
//! @author Dirk Clasen
//! @date 04.08.2005
class cException {
 private:
  std::string m_Message;  ///< error message to be written

 public:
  //! @brief default constructor
  //! @param Message errormessage
  //! @param File name of the file that reported the error, i.e. __FILE__
  //! @param Line row of 'File' where the error occured
  cException(const std::string& Message, const char* File, int Line);
  virtual ~cException();

  //! @brief return the errormessage stored in m_Message
  //! @return the errormessage
  inline const std::string& what(void) const { return m_Message; }
};

#endif
