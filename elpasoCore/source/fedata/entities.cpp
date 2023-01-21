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

#include "entities.h"


void dumpObjectCounters(FILE *os)
{
#ifndef ELPASO_TEST
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  if(tCounter<cNode>::howMany())                               PetscFPrintf(PETSC_COMM_WORLD, os, "  Nodes .......................................: %7d\n", tCounter<cNode>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  if(tCounter<cMaterial>::howMany())                           PetscFPrintf(PETSC_COMM_WORLD, os, "  Materials ...................................: %7d\n", tCounter<cMaterial>::howMany());
  if(tCounter<cMaterialFluid>::howMany())                      PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- fluid ...................................: %7d\n", tCounter<cMaterialFluid>::howMany());
  if(tCounter<cMaterialFluidLin>::howMany())                   PetscFPrintf(PETSC_COMM_WORLD, os, "  |   +-- linear ..............................: %7d\n", tCounter<cMaterialFluidLin>::howMany());
  //if(tCounter<cMaterialFluidIdeal>::howMany())                 PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- ideal fluid .....................: %7d\n", tCounter<cMaterialFluidIdeal>::howMany());
  //if(tCounter<cMaterialFluidLoss>::howMany())                  PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- lossy fluid .....................: %7d\n", tCounter<cMaterialFluidLoss>::howMany());
  //if(tCounter<cMaterialFluidEquiv>::howMany())                 PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- porous equivalent ...............: %7d\n", tCounter<cMaterialFluidEquiv>::howMany());
  //if(tCounter<cMaterialFluidEquivDirect>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- porous equivalent direct.........: %7d\n", tCounter<cMaterialFluidEquivDirect>::howMany());
  // if(tCounter<cMaterialFluidCloaking>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- cloaking fluid ..................: %7d\n", tCounter<cMaterialFluidCloaking>::howMany());
  // if(tCounter<cMaterialPoro>::howMany())                       PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- poroelastic .............................: %7d\n", tCounter<cMaterialPoro>::howMany());
  // if(tCounter<cMaterialPoro>::howMany())                       PetscFPrintf(PETSC_COMM_WORLD, os, "  |   +-- linear ..............................: %7d\n", tCounter<cMaterialPoro>::howMany());
  // if(tCounter<cMaterialStructurePoro>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- poroelastic .....................: %7d\n", tCounter<cMaterialStructurePoro>::howMany());
  // if(tCounter<cMaterialStructurePoroLinear>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "  |       |   +-- linear isotropic ............: %7d\n", tCounter<cMaterialStructurePoroLinear>::howMany());
  // if(tCounter<cMaterialStructurePoroFoam>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "  |       |   +-- foam like ...................: %7d\n", tCounter<cMaterialStructurePoroFoam>::howMany());
  // if(tCounter<cMaterialStructurePoro3d>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- poroelastic 3d ..................: %7d\n", tCounter<cMaterialStructurePoro3d>::howMany());
  // if(tCounter<cMaterialStructurePoro3dtransisotrop>::howMany())PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- poroelastic 3d transv. isotrop ..: %7d\n", tCounter<cMaterialStructurePoro3dtransisotrop>::howMany());
  // if(tCounter<cMaterialStructurePoro3dKappa>::howMany())       PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- poroelastic 3d kappa ............: %7d\n", tCounter<cMaterialStructurePoro3dKappa>::howMany());
  // if(tCounter<cMaterialStructurePoro3dKappaComp>::howMany())   PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- poroelastic 3d kappa comp........: %7d\n", tCounter<cMaterialStructurePoro3dKappaComp>::howMany());
  // if(tCounter<cMaterialStructurePoroPlate>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "  |       +-- poroelastic plate ...............: %7d\n", tCounter<cMaterialStructurePoroPlate>::howMany());
  if(tCounter<cMaterialStructure>::howMany())                  PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure ...............................: %7d\n", tCounter<cMaterialStructure>::howMany());
  if(tCounter<cMaterialMass>::howMany())                       PetscFPrintf(PETSC_COMM_WORLD, os, "      +-- pointmass ...........................: %7d\n", tCounter<cMaterialMass>::howMany());
  if(tCounter<cMaterialStructureIso>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "      +-- isotropic  ..........................: %7d\n", tCounter<cMaterialStructureIso>::howMany());
  if(tCounter<cMaterialStructureIsoLin>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   +-- linear ..........................: %7d\n", tCounter<cMaterialStructureIsoLin>::howMany());
  if(tCounter<cMaterialStructureIsotrop>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- elastic .....................: %7d\n", tCounter<cMaterialStructureIsotrop>::howMany());
  if(tCounter<cMaterialSpring>::howMany())                     PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring ......................: %7d\n", tCounter<cMaterialSpring>::howMany());
  // if(tCounter<cMaterialStructureVisco>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- viscoelastic ................: %7d\n", tCounter<cMaterialStructureVisco>::howMany());
  // if(tCounter<cMaterialStructureViscoCLD_CH>::howMany())       PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- viscoelastic CLD CH .........: %7d\n", tCounter<cMaterialStructureViscoCLD_CH>::howMany());
  // if(tCounter<cMaterialStructureViscoCLD_RKU>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- viscoelastic CLD RKU ........: %7d\n", tCounter<cMaterialStructureViscoCLD_RKU>::howMany());
  // if(tCounter<cMaterialStructureViscoCLD_RKU_Laminate_wa>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- viscoelastic CLD RKU Laminate (w.a.): %7d\n", tCounter<cMaterialStructureViscoCLD_RKU_Laminate_wa>::howMany());
  //if(tCounter<cMaterialStructureViscoConst>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- constant (eta) ..........: %7d\n", tCounter<cMaterialStructureViscoConst>::howMany());
//if(tCounter<cMaterialStructureViscoFrequency>::howMany())    PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- frequency dependent .....: %7d\n", tCounter<cMaterialStructureViscoFrequency>::howMany());
//if(tCounter<cMaterialStructureViscoLFFrequency>::howMany())    PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- frequency dependent LF...: %7d\n", tCounter<cMaterialStructureViscoLFFrequency>::howMany());
//if(tCounter<cMaterialStructureViscoLFFrequencyTvar>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- frequency dependent tvar.: %7d\n", tCounter<cMaterialStructureViscoLFFrequencyTvar>::howMany());
//if(tCounter<cMaterialStructureViscoParametric>::howMany())   PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- parametric description ..: %7d\n", tCounter<cMaterialStructureViscoParametric>::howMany());
  // if(tCounter<cMaterialStructureIsoNonLin>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |   +-- nonlinear .......................: %7d\n", tCounter<cMaterialStructureIsoNonLin>::howMany());
  // if(tCounter<cMaterialStructureIsotropNlEl>::howMany())       PetscFPrintf(PETSC_COMM_WORLD, os, "      |       +-- elastic .....................: %7d\n", tCounter<cMaterialStructureIsotropNlEl>::howMany());
  // if(tCounter<cMaterialStructureIsotropNlEl1D>::howMany())     PetscFPrintf(PETSC_COMM_WORLD, os, "      |       |   +-- 1D ......................: %7d\n", tCounter<cMaterialStructureIsotropNlEl1D>::howMany());
  // if(tCounter<cMaterialStructureElPla>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |       +-- elasticplastic ..............: %7d\n", tCounter<cMaterialStructureElPla>::howMany());
  // if(tCounter<cMaterialCapModel>::howMany())                   PetscFPrintf(PETSC_COMM_WORLD, os, "      |           +-- Cap Model ...............: %7d\n", tCounter<cMaterialCapModel>::howMany());
  // if(tCounter<cMaterialDruckerPrager>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "      |           +-- Drucker Prager ..........: %7d\n", tCounter<cMaterialDruckerPrager>::howMany());
  // if(tCounter<cMaterialMohrCoulomb>::howMany())                PetscFPrintf(PETSC_COMM_WORLD, os, "      |           +-- Mohr Coulomb ............: %7d\n", tCounter<cMaterialMohrCoulomb>::howMany());
  // if(tCounter<cMaterialStructureOrtho>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      +-- orthotropic .........................: %7d\n", tCounter<cMaterialStructureOrtho>::howMany());
  // if(tCounter<cMaterialStructureOrthoLin>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "          +-- linear ..........................: %7d\n", tCounter<cMaterialStructureOrthoLin>::howMany());
  // if(tCounter<cMaterialStructureOrthotrop>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "              +-- elastic .....................: %7d\n", tCounter<cMaterialStructureOrthotrop>::howMany());
  // if(tCounter<cMaterialStructureOrthotropvisco>::howMany())    PetscFPrintf(PETSC_COMM_WORLD, os, "              +-- visco elastic ...............: %7d\n", tCounter<cMaterialStructureOrthotropvisco>::howMany());
  // if(tCounter<cMaterialStructureOrthotropviscoLaminate>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "              +-- visco elastic laminate ......: %7d\n", tCounter<cMaterialStructureOrthotropviscoLaminate>::howMany());
  // if(tCounter<cMaterialStructureOrthotropviscoLaminatePre>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "              +-- visco elastic prestr laminate: %7d\n", tCounter<cMaterialStructureOrthotropviscoLaminatePre>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  if(tCounter<cGroup>::howMany())                              PetscFPrintf(PETSC_COMM_WORLD, os, "  Groups ......................................: %7d\n", tCounter<cGroup>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  if(tCounter<cElement>::howMany())                            PetscFPrintf(PETSC_COMM_WORLD, os, "  Elements ....................................: %7d\n", tCounter<cElement>::howMany());
  if(tCounter<cElementFEM>::howMany())                         PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- FEM .....................................: %7d\n", tCounter<cElementFEM>::howMany());
  if(tCounter<cElementStructure>::howMany())                   PetscFPrintf(PETSC_COMM_WORLD, os, "      +-- structure ...........................: %7d\n", tCounter<cElementStructure>::howMany());
  if(tCounter<cElementStructureMass>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "      |   +-- pointmass .......................: %7d\n", tCounter<cElementStructureMass>::howMany());
  if(tCounter<cElementStructureLinear>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   +-- linear ..........................: %7d\n", tCounter<cElementStructureLinear>::howMany());
  if(tCounter<cElementStructureBeam>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- beam (2 nodes) ..............: %7d\n", tCounter<cElementStructureBeam>::howMany());
  if(tCounter<cElementStructureBeam3D>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- beam3d (2 nodes) ............: %7d\n", tCounter<cElementStructureBeam3D>::howMany());
  if(tCounter<cElementStructureBeam10>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- beam3d (10 dofs) ........: %7d\n", tCounter<cElementStructureBeam10>::howMany());
  if(tCounter<cElementStructureBeam12>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- beam3d (12 dofs) ........: %7d\n", tCounter<cElementStructureBeam12>::howMany());
  if(tCounter<cElementStructureSpring>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring (2 nodes) ............: %7d\n", tCounter<cElementStructureSpring>::howMany());
  if(tCounter<cElementStructureSpringz>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring z (2 nodes) ..........: %7d\n", tCounter<cElementStructureSpringz>::howMany());
  if(tCounter<cElementStructureSpringBC>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring bc (1 node) ..........: %7d\n", tCounter<cElementStructureSpringBC>::howMany());
  if(tCounter<cElementStructureSpringBCx>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring bc x (1 node) ........: %7d\n", tCounter<cElementStructureSpringBCx>::howMany());
  if(tCounter<cElementStructureSpringBCy>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring bc y (1 node) ........: %7d\n", tCounter<cElementStructureSpringBCy>::howMany());
  if(tCounter<cElementStructureSpringBCz>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring bc z (1 node) ........: %7d\n", tCounter<cElementStructureSpringBCz>::howMany());
  if(tCounter<cElementStructureSpringBCrx>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring bc rx (1 node) .......: %7d\n", tCounter<cElementStructureSpringBCrx>::howMany());
  if(tCounter<cElementStructureSpringBCry>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring bc ry (1 node) .......: %7d\n", tCounter<cElementStructureSpringBCry>::howMany());
  if(tCounter<cElementStructureSpringBCrz>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- spring bc rz (1 node) .......: %7d\n", tCounter<cElementStructureSpringBCrz>::howMany());
  // if(tCounter<cElementStructureDisc>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- membrane ....................: %7d\n", tCounter<cElementStructureDisc>::howMany());
  // if(tCounter<cElementStructureDisc3>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- membrane (3 nodes) ......: %7d\n", tCounter<cElementStructureDisc3>::howMany());
  // if(tCounter<cElementStructureDisc4>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- membrane (4 nodes) ......: %7d\n", tCounter<cElementStructureDisc4>::howMany());
  // if(tCounter<cElementStructureDisc9>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- membrane (9 nodes) ......: %7d\n", tCounter<cElementStructureDisc9>::howMany());
  // if(tCounter<cElementStructureDiscDrill>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- membrane w. phi_z ...........: %7d\n", tCounter<cElementStructureDiscDrill>::howMany());
  // if(tCounter<cElementStructureDiscDrill4>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- membrane (4 nodes) ......: %7d\n", tCounter<cElementStructureDiscDrill4>::howMany());
  // if(tCounter<cElementStructureMindlin>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- Mindlin .....................: %7d\n", tCounter<cElementStructureMindlin>::howMany());
  // if(tCounter<cElementStructureMindlinDSG3>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- DSG (3 nodes) ...........: %7d\n", tCounter<cElementStructureMindlinDSG3>::howMany());
  // if(tCounter<cElementStructureMindlinDSG4>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- DSG (4 nodes) ...........: %7d\n", tCounter<cElementStructureMindlinDSG4>::howMany());
  // if(tCounter<cElementStructureMindlinDSG9>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- DSG (9 nodes) ...........: %7d\n", tCounter<cElementStructureMindlinDSG9>::howMany());
  // if(tCounter<cElementStructureMindlinDSG9pre>::howMany())     PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- DSG prestressed (9 nodes): %7d\n", tCounter<cElementStructureMindlinDSG9pre>::howMany());
  // if(tCounter<cElementStructureKienzler>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- Kienzler ....................: %7d\n", tCounter<cElementStructureKienzler>::howMany());
  if(tCounter<cElementStructureKirchhoff>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- Kirchhoff ...................: %7d\n", tCounter<cElementStructureKirchhoff>::howMany());
  if(tCounter<cElementStructureKirchhoff4>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- 4 nodes .................: %7d\n", tCounter<cElementStructureKirchhoff4>::howMany());
  // if(tCounter<cElementStructurePlaneShell>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- plane shell element .........: %7d\n", tCounter<cElementStructurePlaneShell>::howMany());
  // if(tCounter<cElementStructurePlaneShell3>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- 3 nodes .................: %7d\n", tCounter<cElementStructurePlaneShell3>::howMany());
  // if(tCounter<cElementStructurePlaneShell4>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- 4 nodes .................: %7d\n", tCounter<cElementStructurePlaneShell4>::howMany());
  // if(tCounter<cElementStructurePlaneShell9>::howMany())        PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- 9 nodes .................: %7d\n", tCounter<cElementStructurePlaneShell9>::howMany());
  // if(tCounter<cElementStructurePlaneShell9pre>::howMany())     PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- 9 nodes prestressed......: %7d\n", tCounter<cElementStructurePlaneShell9pre>::howMany());
  // if(tCounter<cElementStructurePlaneShellDrill>::howMany())    PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- pl. shell elem. (drill.) ....: %7d\n", tCounter<cElementStructurePlaneShellDrill>::howMany());
  // if(tCounter<cElementStructurePlaneShellDrill4>::howMany())   PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- 4 nodes .................: %7d\n", tCounter<cElementStructurePlaneShellDrill4>::howMany());
  // if(tCounter<cElementStructureTetra>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- tetra .......................: %7d\n", tCounter<cElementStructureTetra>::howMany());
  // if(tCounter<cElementStructureTetra4>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- tetra (4 nodes) .........: %7d\n", tCounter<cElementStructureTetra4>::howMany());
  // if(tCounter<cElementStructureTetra4L>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- tetraL(4 nodes) .........: %7d\n", tCounter<cElementStructureTetra4L>::howMany());
  // if(tCounter<cElementStructureTetra10>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- tetra (10 nodes).........: %7d\n", tCounter<cElementStructureTetra10>::howMany());
  // if(tCounter<cElementStructureTetra10L>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- tetraL(10 nodes).........: %7d\n", tCounter<cElementStructureTetra10L>::howMany());
  // if(tCounter<cElementStructureTetra16>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- tetra (16 nodes).........: %7d\n", tCounter<cElementStructureTetra16>::howMany());
  // if(tCounter<cElementStructureTetra16L>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- tetraL(16 nodes).........: %7d\n", tCounter<cElementStructureTetra16L>::howMany());
  if(tCounter<cElementStructureBrick>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- brick .......................: %7d\n", tCounter<cElementStructureBrick>::howMany());
  if(tCounter<cElementStructureBrick8>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- brick ( 8 nodes) ........: %7d\n", tCounter<cElementStructureBrick8>::howMany());
  //if(tCounter<cElementStructureBrick20>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- brick (20 nodes) ........: %7d\n", tCounter<cElementStructureBrick20>::howMany());
  //if(tCounter<cElementStructureBrick27>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- brick (27 nodes) ........: %7d\n", tCounter<cElementStructureBrick27>::howMany());
  // if(tCounter<cElementStructureNonLinear>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |   +-- nonlinear .......................: %7d\n", tCounter<cElementStructureNonLinear>::howMany());
  // if(tCounter<cElementStructureBeam3DNL>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- beam3d (2 nodes) ............: %7d\n", tCounter<cElementStructureBeam3DNL>::howMany());
  // if(tCounter<cElementStructureBeam10NL>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- beam3d (10 dofs) ........: %7d\n", tCounter<cElementStructureBeam10NL>::howMany());
  // if(tCounter<cElementStructureCable>::howMany())              PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- cable .......................: %7d\n", tCounter<cElementStructureCable>::howMany());
  // if(tCounter<cElementStructureCable2D>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- cable2d (2 nodes) .......: %7d\n", tCounter<cElementStructureCable2D>::howMany());
  // if(tCounter<cElementStructureCable3D>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- cable3d (2 nodes) .......: %7d\n", tCounter<cElementStructureCable3D>::howMany());
  // if(tCounter<cElementStructurePoro>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "      |   +-- poro ............................: %7d\n", tCounter<cElementStructurePoro>::howMany());
  // if(tCounter<cElementStructurePoroPlate>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- plate .......................: %7d\n", tCounter<cElementStructurePoroPlate>::howMany());
  // if(tCounter<cElementStructurePoroPlateQ2P1>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- Q2P1 (4 nodes) ..........: %7d\n", tCounter<cElementStructurePoroPlateQ2P1>::howMany());
  // if(tCounter<cElementStructurePoroPlateQ3P1>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- Q3P1 (9 nodes) ..........: %7d\n", tCounter<cElementStructurePoroPlateQ3P1>::howMany());
  // if(tCounter<cElementStructurePoroPlateQ2P2>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- Q2P2 (4 nodes) ..........: %7d\n", tCounter<cElementStructurePoroPlateQ2P2>::howMany());
  // if(tCounter<cElementStructurePoroPlateQ2P3>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- Q2P3 (4 nodes) ..........: %7d\n", tCounter<cElementStructurePoroPlateQ2P3>::howMany());
  // if(tCounter<cElementStructurePoroPlateQ3P3>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- Q3P3 (9 nodes) ..........: %7d\n", tCounter<cElementStructurePoroPlateQ3P3>::howMany());
  // if(tCounter<cElementStructurePoroPlateKienzler>::howMany())  PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- PoroPlateKienzler .......: %7d\n", tCounter<cElementStructurePoroPlateKienzler>::howMany());
  // if(tCounter<cElementStructurePoroPlateKienzler4>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   |   +-- PoroPlateKienzler4 ..: %7d\n", tCounter<cElementStructurePoroPlateKienzler4>::howMany());
  // if(tCounter<cElementStructurePoroDisc>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- disc ........................: %7d\n", tCounter<cElementStructurePoroDisc>::howMany());
  // if(tCounter<cElementStructurePoroDiscQ2P1>::howMany())       PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   |   +-- Q2P1 (4 nodes) ..........: %7d\n", tCounter<cElementStructurePoroDiscQ2P1>::howMany());
  // if(tCounter<cElementStructurePoroDiscKienzler>::howMany())  PetscFPrintf(PETSC_COMM_WORLD, os,  "      |   |   |   +-- PoroDiscKienzler ........: %7d\n", tCounter<cElementStructurePoroDiscKienzler>::howMany());
  // if(tCounter<cElementStructurePoroDiscKienzler4>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os,  "      |   |   |   |   +-- PoroDiscKienzler4 ...: %7d\n", tCounter<cElementStructurePoroDiscKienzler4>::howMany());
  // if(tCounter<cElementStructurePoroPlaneShell>::howMany())     PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |   +-- plane shell .................: %7d\n", tCounter<cElementStructurePoroPlaneShell>::howMany());
  // if(tCounter<cElementStructurePoroPlaneShellQ2P1>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |       +-- Q2P1 (4 nodes) ..........: %7d\n", tCounter<cElementStructurePoroPlaneShellQ2P1>::howMany());
  // if(tCounter<cElementStructurePoroPlaneShellKienzler>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |    +-- PoroShellKienzler ..........: %7d\n", tCounter<cElementStructurePoroPlaneShellKienzler>::howMany());
  // if(tCounter<cElementStructurePoroPlaneShellKienzler4>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "      |   |    |   +-- PoroShellKienzler4 .....: %7d\n", tCounter<cElementStructurePoroPlaneShellKienzler4>::howMany());
  // if(tCounter<cElementStructurePoro3d>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "      |   +-- poro 3d .........................: %7d\n", tCounter<cElementStructurePoro3d>::howMany());
  // if(tCounter<cElementStructurePoro3dUP>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |       +-- u-p formulation .............: %7d\n", tCounter<cElementStructurePoro3dUP>::howMany());
  // if(tCounter<cElementStructurePoro3dUP8>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |       |   +-- Poro3dUP8   (8 nodes) ...: %7d\n", tCounter<cElementStructurePoro3dUP8>::howMany());
  // if(tCounter<cElementStructurePoro3dUP27>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "      |       |   +-- Poro3dUP27 (27 nodes) ...: %7d\n", tCounter<cElementStructurePoro3dUP27>::howMany());
  // if(tCounter<cElementStructurePoro3dUU>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "      |       +-- u-U formulation .............: %7d\n", tCounter<cElementStructurePoro3dUU>::howMany());
  // if(tCounter<cElementStructurePoro3dUU8>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "      |           +-- Poro3dUU8   (8 nodes) ...: %7d\n", tCounter<cElementStructurePoro3dUU8>::howMany());
  if(tCounter<cElementFluid>::howMany())                       PetscFPrintf(PETSC_COMM_WORLD, os, "      +-- fluid ...............................: %7d\n", tCounter<cElementFluid>::howMany());
  if(tCounter<cElementFluid3d>::howMany())                     PetscFPrintf(PETSC_COMM_WORLD, os, "          +-- 3D ..............................: %7d\n", tCounter<cElementFluid3d>::howMany());
  // if(tCounter<cElementFluid4>::howMany())                      PetscFPrintf(PETSC_COMM_WORLD, os, "          |   +-- fluid ( 4 nodes) ............: %7d\n", tCounter<cElementFluid4>::howMany());
  if(tCounter<cElementFluid8>::howMany())                      PetscFPrintf(PETSC_COMM_WORLD, os, "          |   +-- fluid ( 8 nodes) ............: %7d\n", tCounter<cElementFluid8>::howMany());
  // if(tCounter<cElementFluid27>::howMany())                     PetscFPrintf(PETSC_COMM_WORLD, os, "          |   +-- fluid (27 nodes) ............: %7d\n", tCounter<cElementFluid27>::howMany());
  // if(tCounter<cElementFluid2d>::howMany())                     PetscFPrintf(PETSC_COMM_WORLD, os, "          +-- 2D ..............................: %7d\n", tCounter<cElementFluid2d>::howMany());
  // if(tCounter<cElementFluid2d4>::howMany())                    PetscFPrintf(PETSC_COMM_WORLD, os, "              +-- fluid ( 4 nodes) ............: %7d\n", tCounter<cElementFluid2d4>::howMany());
  // if(tCounter<cElementFluid2d9>::howMany())                    PetscFPrintf(PETSC_COMM_WORLD, os, "              +-- fluid ( 9 nodes) ............: %7d\n", tCounter<cElementFluid2d9>::howMany());
  // if(tCounter<cElementFF>::howMany())                          PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- fluid flow (auxiliary elems) ............: %7d\n", tCounter<cElementFF>::howMany());
  // if(tCounter<cElementFF4>::howMany())                         PetscFPrintf(PETSC_COMM_WORLD, os, "          +-- FF4 ( 4 nodes) ..................: %7d\n", tCounter<cElementFF4>::howMany());
  // if(tCounter<cElementFF9>::howMany())                         PetscFPrintf(PETSC_COMM_WORLD, os, "          +-- FF9 ( 9 nodes) ..................: %7d\n", tCounter<cElementFF9>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " ---------------------------------------------------------\n");
  if(tCounter<cElementInterface>::howMany())                   PetscFPrintf(PETSC_COMM_WORLD, os, "  Interface elements ..........................: %7d\n", tCounter<cElementInterface>::howMany());
  // if(tCounter<cElementInterfaceBeam>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- 2d acoustic + Timoshenko beams ..........: %7d\n", tCounter<cElementInterfaceBeam>::howMany());
  // if(tCounter<cElementInterfaceKirchhoff>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- acoustic + Kirchhoff plate ..............: %7d\n", tCounter<cElementInterfaceKirchhoff>::howMany());
  if(tCounter<cElementInterfaceMindlin>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- acoustic + Mindlin plate ................: %7d\n", tCounter<cElementInterfaceMindlin>::howMany());
  // if(tCounter<cElementInterfacePoro>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- acoustic + Poroelastic plate ............: %7d\n", tCounter<cElementInterfacePoro>::howMany());
  // if(tCounter<cElementInterfaceMindlinPoro3dUP>::howMany())    PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- Poro3dUP + Mindlin plate (bounded) ......: %7d\n", tCounter<cElementInterfaceMindlinPoro3dUP>::howMany());
  // if(tCounter<cElementInterfaceMindlinEquiporo>::howMany())    PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- equiporo + Mindlin plate (bounded).......: %7d\n", tCounter<cElementInterfaceMindlinEquiporo>::howMany());
  // if(tCounter<cElementInterfaceFFPoro3dUP>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- fluid flow + Poro3dUP ...................: %7d\n", tCounter<cElementInterfaceFFPoro3dUP>::howMany());
  // if(tCounter<cElementInterfaceFFShell>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- fluid flow + Shell ......................: %7d\n", tCounter<cElementInterfaceFFShell>::howMany());
  // if(tCounter<cElementInterfacePlaneShellPoro2dUP>::howMany()) PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- PoroShellKienzler4 + plane shell ........: %7d\n", tCounter<cElementInterfacePlaneShellPoro2dUP>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " ---------------------------------------------------------\n");
  if(tCounter<cElementCoupling>::howMany())                    PetscFPrintf(PETSC_COMM_WORLD, os, "  FEM-BEM-coupling ............................: %7d\n", tCounter<cElementCoupling>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  if(tCounter<cBoundaryCondition>::howMany())                  PetscFPrintf(PETSC_COMM_WORLD, os, "  Nodal boundary conditions ...................: %7d\n", tCounter<cBoundaryCondition>::howMany());
  if(tCounter<cBoundaryConditionStructure>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure ...............................: %7d\n", tCounter<cBoundaryConditionStructure>::howMany());
  if(tCounter<cBoundaryConditionFluid>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- fluid ...................................: %7d\n", tCounter<cBoundaryConditionFluid>::howMany());
  // if(tCounter<cBoundaryConditionFluidOblique>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "      +-- oblique incidence ...................: %7d\n", tCounter<cBoundaryConditionFluidOblique>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  // if(tCounter<cConstraint>::howMany())                         PetscFPrintf(PETSC_COMM_WORLD, os, "  Constraints .................................: %7d\n", tCounter<cConstraint>::howMany());
  // if(tCounter<cConstraintMultiPoint>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- multi point constraints..................: %7d\n", tCounter<cConstraintMultiPoint>::howMany());
  // if(tCounter<cConstraintMeshTie>::howMany())                  PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- mesh tie constraints.....................: %7d\n", tCounter<cConstraintMeshTie>::howMany());
  // PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  if(tCounter<cNodalForce>::howMany())                         PetscFPrintf(PETSC_COMM_WORLD, os, "  Nodal forces ................................: %7d\n", tCounter<cNodalForce>::howMany());
  if(tCounter<cNodalForceStructure>::howMany())                PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure ...............................: %7d\n", tCounter<cNodalForceStructure>::howMany());
  // if(tCounter<cNodalForceStructureTime>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure (time dependent) ..............: %7d\n", tCounter<cNodalForceStructureTime>::howMany());
  if(tCounter<cNodalForceFluid>::howMany())                    PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- fluid ...................................: %7d\n", tCounter<cNodalForceFluid>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  if(tCounter<cNodalMoment>::howMany())                        PetscFPrintf(PETSC_COMM_WORLD, os, "  Nodal moments ...............................: %7d\n", tCounter<cNodalMoment>::howMany());
  if(tCounter<cNodalMomentStructure>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure ...............................: %7d\n", tCounter<cNodalMomentStructure>::howMany());
  // if(tCounter<cNodalMomentStructureTime>::howMany())           PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure (time dependent) ..............: %7d\n", tCounter<cNodalMomentStructureTime>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
//   if(tCounter<cNodalPressure>::howMany())                      PetscFPrintf(PETSC_COMM_WORLD, os, "  Nodal pressures .............................: %7d\n", tCounter<cNodalPressure>::howMany());
//   if(tCounter<cNodalPressureStructure>::howMany())             PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure ...............................: %7d\n", tCounter<cNodalPressureStructure>::howMany());
// #ifdef HAVE_FFTW
//   if(tCounter<cNodalPressureStructureFreq>::howMany())         PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure (frequency dependent) .........: %7d\n", tCounter<cNodalPressureStructureFreq>::howMany());
// #endif
//   PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
//   if(tCounter<cNodalValues>::howMany())                        PetscFPrintf(PETSC_COMM_WORLD, os, "  Nodal values ................................: %7d\n", tCounter<cNodalValues>::howMany());
// #ifdef HAVE_FFTW
//   if(tCounter<cNodalValuesFF>::howMany())                      PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- values ff (fluid flow)...................: %7d\n", tCounter<cNodalValuesFF>::howMany());
//   if(tCounter<cNodalValuesVelocity>::howMany())                 PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- values velocity .........................: %7d\n", tCounter<cNodalValuesVelocity>::howMany());

// #endif
  // PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");

  if(tCounter<cElementLoad>::howMany())                        PetscFPrintf(PETSC_COMM_WORLD, os, "  Elementloads ................................: %7d\n", tCounter<cElementLoad>::howMany());
  if(tCounter<cElementLoadStructure>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- structure................................: %7d\n", tCounter<cElementLoadStructure>::howMany());
  if(tCounter<cElementLoadStructureConst>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "  |   +-- constant ............................: %7d\n", tCounter<cElementLoadStructureConst>::howMany());
  // if(tCounter<cElementLoadStructureFrq>::howMany())            PetscFPrintf(PETSC_COMM_WORLD, os, "  |   +-- frequency-dependent..................: %7d\n", tCounter<cElementLoadStructureFrq>::howMany());
  // if(tCounter<cElementLoadStructureDistributed>::howMany())    PetscFPrintf(PETSC_COMM_WORLD, os, "  |   +-- varying .............................: %7d\n", tCounter<cElementLoadStructureDistributed>::howMany());
  // if(tCounter<cElementLoadVn>::howMany())                      PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- normal sound velocity ...................: %7d\n", tCounter<cElementLoadVn>::howMany());
  // if(tCounter<cElementLoadImpedance>::howMany())               PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- impedance ...............................: %7d\n", tCounter<cElementLoadImpedance>::howMany());
  // if(tCounter<cElementLoadImpedanceConst>::howMany())          PetscFPrintf(PETSC_COMM_WORLD, os, "  |   +-- const ...............................: %7d\n", tCounter<cElementLoadImpedanceConst>::howMany());
  // if(tCounter<cElementLoadImpedanceFrequency>::howMany())      PetscFPrintf(PETSC_COMM_WORLD, os, "  |   +-- frequency ...........................: %7d\n", tCounter<cElementLoadImpedanceFrequency>::howMany());
  // if(tCounter<cElementLoadSBFEM>::howMany())                   PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- SBFEM ...................................: %7d\n", tCounter<cElementLoadSBFEM>::howMany());
  PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
  // if(tCounter<cInterface>::howMany())                          PetscFPrintf(PETSC_COMM_WORLD, os, "  Inferfaces ..................................: %7d\n", tCounter<cInterface>::howMany());
  // if(tCounter<cSwebem>::howMany())                             PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- SWEBEM ..................................: %7d\n", tCounter<cSwebem>::howMany());
  // if(tCounter<cBem>::howMany())                                PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- tBEM ....................................: %7d\n", tCounter<cBem>::howMany());
  // if(tCounter<cSimilar>::howMany())                            PetscFPrintf(PETSC_COMM_WORLD, os, "  +-- ScaBo ...................................: %7d\n", tCounter<cSimilar>::howMany());
  // PetscFPrintf(PETSC_COMM_WORLD, os, " --------------------------------------------------------- \n");
#endif // ELPASO_TEST
}
