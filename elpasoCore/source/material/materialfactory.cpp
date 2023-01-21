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

#include "materialfactory.h"
#include "material.h"
#include "../misc/parser/femparserinterface.h"

cMaterialFactory::cMaterialFactory()
{
    // empty
}

cMaterialFactory::~cMaterialFactory()
{
    // empty
}

cMaterial* cMaterialFactory::createMaterial(ParserMaterialData _matData)
{
    std::stringstream   data;
    cMaterial*          pMaterial = 0;
    // -------------------------------------------------------------------------
    //   create a proper new object according to the attribute Typ
    // -------------------------------------------------------------------------
    if (_matData.matType == "STR_LIN_ELA_ISO_DIR")
    {
        pMaterial = new cMaterialStructureIsotrop;
        data << _matData.matParameters[(eMaterialTags)Id] << " ";
        data << _matData.matParameters[(eMaterialTags)E] << " ";
        data << _matData.matParameters[(eMaterialTags)nu] << " ";
        data << _matData.matParameters[(eMaterialTags)A] << " ";
        data << _matData.matParameters[(eMaterialTags)Ix] << " ";
        data << _matData.matParameters[(eMaterialTags)Iy] << " ";
        data << _matData.matParameters[(eMaterialTags)Iz] << " ";
        data << _matData.matParameters[(eMaterialTags)rho] << " ";
        data << _matData.matParameters[(eMaterialTags)t] << " ";
        data << _matData.matParameters[(eMaterialTags)Fi] << " ";
        //data >> *pMaterial;
    }
    else if (_matData.matType == "STR_LIN_SPR_ORT_DIR")
    {
        pMaterial = new cMaterialSpring;
        data << _matData.matParameters[(eMaterialTags)Id] << " ";
        data << _matData.matParameters[(eMaterialTags)Cx] << " ";
        data << _matData.matParameters[(eMaterialTags)Cy] << " ";
        data << _matData.matParameters[(eMaterialTags)Cz] << " ";
        data << _matData.matParameters[(eMaterialTags)Crx] << " ";
        data << _matData.matParameters[(eMaterialTags)Cry] << " ";
        data << _matData.matParameters[(eMaterialTags)Crz] << " ";
        data << _matData.matParameters[(eMaterialTags)eta] << " ";
        data << _matData.matParameters[(eMaterialTags)mtyp] << " ";
        //data >> *pMaterial;
    }
    else if (_matData.matType == "pointmass")
    {
        pMaterial = new cMaterialMass;
        data << _matData.matParameters[(eMaterialTags)Id] << " ";
        data << _matData.matParameters[(eMaterialTags)M] << " ";
        //data >> *pMaterial;
    }
    else if (_matData.matType == "AF_LIN_UAF_ISO_DIR")
    {
        pMaterial = new cMaterialFluidIdeal;
        data << _matData.matParameters[(eMaterialTags)Id] << " ";
        data << _matData.matParameters[(eMaterialTags)c] << " ";
        data << _matData.matParameters[(eMaterialTags)rho] << " ";
        data << _matData.matParameters[(eMaterialTags)t] << " ";
        //data >> *pMaterial;
    }
    else if (_matData.matType == "STR_LIN_VIS_ISO_DIR")
    {
        pMaterial = new cMaterialStructureVisco;
        data << _matData.matParameters[(eMaterialTags)Id] << " ";
        data << _matData.matParameters[(eMaterialTags)mtyp].c_str() << " ";
        data << _matData.matParameters[(eMaterialTags)E] << " ";
        data << _matData.matParameters[(eMaterialTags)eta] << " ";
        data << _matData.matParameters[(eMaterialTags)nu] << " ";
        data << _matData.matParameters[(eMaterialTags)A] << " ";
        data << _matData.matParameters[(eMaterialTags)Ix] << " ";
        data << _matData.matParameters[(eMaterialTags)Iy] << " ";
        data << _matData.matParameters[(eMaterialTags)Iz] << " ";
        data << _matData.matParameters[(eMaterialTags)rho] << " ";
        data << _matData.matParameters[(eMaterialTags)t] << " ";
        //data >> *pMaterial;
    }
    else
    {
        return pMaterial;
    }

    // -------------------------------------------------------------------------
    //   now insert the data into the material object
    // -------------------------------------------------------------------------
    data >> *pMaterial;
    pMaterial->setIdentifier(_matData.matName);
    return pMaterial;
}