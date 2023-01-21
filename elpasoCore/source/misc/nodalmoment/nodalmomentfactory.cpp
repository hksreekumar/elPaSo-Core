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

#include "nodalmomentfactory.h"
#include "../parser/femparserinterface.h"
#include "nodalmomentstructure.h"

cNodalMomentFactory::cNodalMomentFactory()
{
    // empty
}

cNodalMomentFactory::~cNodalMomentFactory()
{
    // empty
}

cNodalMoment* cNodalMomentFactory::createNodalMoment(ParserNodeMomentsData _data)
{
    cNodalMoment*           ptr = 0;
    std::stringstream       data_stream;

    
    if (_data.nmomentType == "point_moment")
    {
        ptr = new cNodalMomentStructure;
        data_stream << _data.nmomentId << " " << _data.nmomentMx << " " << _data.nmomentMy << " " << _data.nmomentMz << " ";
        data_stream >> *ptr;
    }

    return ptr;
}