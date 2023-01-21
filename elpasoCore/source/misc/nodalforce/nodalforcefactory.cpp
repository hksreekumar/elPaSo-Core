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

#include "nodalforcefactory.h"
#include "nodalforcestructure.h"
#include "nodalforcefluid.h"
#include "../parser/femparserinterface.h"

cNodalForceFactory::cNodalForceFactory()
{
    // empty
}

cNodalForceFactory::~cNodalForceFactory()
{
    // empty
}

cNodalForce* cNodalForceFactory::createNodalForce(ParserNodeLoadsData _data)
{
    cNodalForce*        ptr = 0;
    std::stringstream   data_stream;

    if (_data.nloadType == "point_force")
    {
        ptr = new cNodalForceStructure;
        data_stream << _data.nloadId << " " << _data.nloadFx << " " << _data.nloadFy << " " << _data.nloadFz << " ";
        data_stream >> *ptr;
    }
    else if (_data.nloadType == "fluid")
    {
        ptr = new cNodalForceFluid;
        data_stream << _data.nloadId << " " << _data.nloadPf << " ";
        data_stream >> *ptr;
    }

    return ptr;
}