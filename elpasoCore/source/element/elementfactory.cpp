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

#include "elementfactory.h"
#include "elementfem.h"
#include "../fedata/entities.h"

cElementFactory::cElementFactory()
{
    // empty
}

cElementFactory::~cElementFactory()
{
    // empty
}

cElementFEM* cElementFactory::createElement(std::string _elementType)
{
    cElementFEM* pElement = 0;

    //****************************************************************************************
    //                             StructureLinear
    //****************************************************************************************
    // <Beam> element (Euler-Bernoulli theory)
    if (_elementType == "Beam")
        pElement = new cElementStructureBeam(Bernoulli);
    // <Beam> element (Euler-Bernoulli theory)
    else if (_elementType == "BeamBernoulli")
        pElement = new cElementStructureBeam(Bernoulli);
    // <BeamBernoulli10> element (Euler-Bernoulli theory) 10 dofs
    else if (_elementType == "BeamBernoulli10")
        pElement = new cElementStructureBeam10(Bernoulli);
    // <BeamBernoulli3d12> element (Euler-Bernoulli theory) 12 dofs
    else if (_elementType == "BeamBernoulli12")
        pElement = new cElementStructureBeam12(Bernoulli);
    // <BeamTimoshenko> element (Timoshenko theory)
    else if (_elementType == "BeamTimoshenko")
        pElement = new cElementStructureBeam(Timoshenko);
    // <BeamTimoshenko10> element (Timoshenko theory) 10 dofs
    else if (_elementType == "BeamTimoshenko10")
        pElement = new cElementStructureBeam10(Timoshenko);
    // <BeamTimoshenko3d12> element (Timoshenko theory) 12 dofs
    else if (_elementType == "BeamTimoshenko12")
        pElement = new cElementStructureBeam12(Timoshenko);
    // <Pointmass> element
    else if (_elementType == "Pointmass")
        pElement = new cElementStructureMass;
    // <Spring> element
    else if (_elementType == "Spring")
        pElement = new cElementStructureSpring;
    // <Springz> element
    else if (_elementType == "Springz")
        pElement = new cElementStructureSpringz;
    // <SpringBC> element
    else if (_elementType == "SpringBC")
        pElement = new cElementStructureSpringBC;
    // <SpringBCx> element
    else if (_elementType == "SpringBCx")
        pElement = new cElementStructureSpringBCx;
    // <SpringBCy> element
    else if (_elementType == "SpringBCy")
        pElement = new cElementStructureSpringBCy;
    // <SpringBCz> element
    else if (_elementType == "SpringBCz")
        pElement = new cElementStructureSpringBCz;
    // <SpringBCrx> element
    else if (_elementType == "SpringBCrx")
        pElement = new cElementStructureSpringBCrx;
    // <SpringBCry> element
    else if (_elementType == "SpringBCry")
        pElement = new cElementStructureSpringBCry;
    // <SpringBCrz> element
    else if (_elementType == "SpringBCrz")
        pElement = new cElementStructureSpringBCrz;
    // <Brick8> element
    else if (_elementType == "Brick8")
        pElement = new cElementStructureBrick8;
    // Mindlin plate element with 4 nodes (<DSG4>)
    else if (_elementType == "DSG4")
        pElement = new cElementStructureMindlinDSG4;
    // Kirchhoff plate element with 4 nodes (<Kirch4>)
    else if (_elementType == "Kirch4")
        pElement = new cElementStructureKirchhoff4;
    //****************************************************************************************
    //                             FluidLinear
    //****************************************************************************************
    // <Fluid8> element (linear hexahedron)
    else if (_elementType == "Fluid8")
        pElement = new cElementFluid8;

    return pElement;
}