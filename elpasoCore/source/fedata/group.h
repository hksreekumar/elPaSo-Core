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

#ifndef INFAM_GROUP_H
#define INFAM_GROUP_H

#include "../misc/counter.h"
#include "../misc/id.h"

//! @brief stores elements and the material assigned to them
//! @author Dirk Clasen
//! @date 15.11.2006
//! stores all data that has been defined in Patran for a specific group
class cGroup : public cId, private tCounter<cGroup> {
 private:
  std::string m_Identifier;   ///< identifier of group
  PetscInt m_MaterialId;      ///< id of material assigned to group
  std::string m_Orientation;  ///< orientation of the material in all elements
                              ///< in that group (local coord, global coord or
                              ///< user-defined by external file)
  std::string
      m_OrientationFilename;  ///< filename giving orientation of the material
                              ///< in all elements in that group
  std::set<PetscInt> m_ElementIds;  ///< ids of elements in group

 public:
  cGroup();
  ~cGroup();

  //! return the number of instances of this object
  static size_t howMany(void) { return tCounter<cGroup>::howMany(); }

  //! assign a new identifier to group
  void setIdentifier(const std::string &ident) { m_Identifier = ident; }

  //! return identifier of group
  std::string getIdentifier(void) const { return m_Identifier; }

  //! read id of material assigned to group
  PetscInt getMaterialId(void) const { return m_MaterialId; }

  //! read orientation of the material assigned to group
  std::string getOrientation(void) const { return m_Orientation; }

  //! read orientation of the material assigned to group
  std::string getOrientationFilename(void) const {
    return m_OrientationFilename;
  }

  //! assign new material to group
  //! @param Id id of material to assign to group
  void setMaterialId(const PetscInt &Id) { m_MaterialId = Id; }

  //! set the orientation for the group
  //! @param Orientation of material in each element of the group
  void setOrientation(const std::string &orientation) {
    m_Orientation = orientation;
  }

  //! set the orientation for the group
  //! @param Orientation of material in each element of the group
  void setOrientationFilename(const std::string &filename) {
    m_OrientationFilename = filename;
  }

  //! insert new element into group
  //! @param Id id of element to insert
  void insertElement(const PetscInt &Id) { m_ElementIds.insert(Id); }

  //! delete an element
  void eraseElement(const PetscInt &Id) { m_ElementIds.erase(Id); }

  std::set<PetscInt>::iterator getFirstElementId(void) {
    return m_ElementIds.begin();
  }
  std::set<PetscInt>::iterator getLastElementId(void) {
    return m_ElementIds.end();
  }

  PetscInt getNumberOfElementsInGroup(void) const {
    return (PetscInt)m_ElementIds.size();
  }
};

#endif
