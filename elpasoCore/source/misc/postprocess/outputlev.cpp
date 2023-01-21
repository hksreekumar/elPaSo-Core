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

#include "outputlev.h"
//#include "../element/fluid/elementfluid3d.h"

cOutputLEV::cOutputLEV()
{  

  this->m_GlobalNumRowsPerRoom[0] = 0;
  this->m_GlobalNumRowsPerRoom[1] = 0;

  this->m_LevelFileName.clear();
}

cOutputLEV::cOutputLEV(cOutputLEV &other):
cOutputfile(other)
{
  this->m_GlobalNumRowsPerRoom[0] = other.m_GlobalNumRowsPerRoom[0];
  this->m_GlobalNumRowsPerRoom[1] = other.m_GlobalNumRowsPerRoom[1];

  this->m_LevelFileName.clear();
  //this->m_LevelFileName	= std::string(other.m_LevelFileName);
}

cOutputLEV::~cOutputLEV()
{
  if (wantOutput() > 0) {
    m_LevelFile.close();
  }
}


void cOutputLEV::writeResultsToFile(cMesh &MyMesh, const int &NrStep, const double &Step)
{
  if (wantOutput() > 0) {

    computePressureLevel(Step, MyMesh);
  }
}


void cOutputLEV::checkForSides(cMesh &MyMesh)
{
  // -------------------------------------------------------------------------
  //   create file name
  // -------------------------------------------------------------------------
  std::string text = MyMesh.getFilename();
  if (text.rfind(cstInputFileSuffix) == (text.length()-cstInputFileSuffix.length()))
    text.erase(text.length()-cstInputFileSuffix.length());
  m_LevelFileName = text + ".lev";


  for (ItMapNodes it = MyMesh.getFirstNode(); it != MyMesh.getLastNode(); it++)
  {
    if (it->second->getRoom() == 1)
      m_GlobalNumRowsPerRoom[0]++;
    else if (it->second->getRoom() == 2)
      m_GlobalNumRowsPerRoom[1]++;
  }

  message("  total number : %6d\n", m_GlobalNumRowsPerRoom[0]+m_GlobalNumRowsPerRoom[1]);
  message("        room A : %6d\n", m_GlobalNumRowsPerRoom[0]);
  message("        room B : %6d\n", m_GlobalNumRowsPerRoom[1]);

  trace(" finished!");
}



void cOutputLEV::parseNdsFile(cMesh &MyMesh)
{
  std::ifstream rein;  // inputstream
  std::string   text;

  short seite   = 1;     // side [1,2]


  // -------------------------------------------------------------------------
  //   opening nds file
  //   used to determine the nodes belonging to the two adjacent rooms
  // -------------------------------------------------------------------------
  text = MyMesh.getFilename();
  if (text.rfind(cstInputFileSuffix) == (text.length()-cstInputFileSuffix.length()))
    text.erase(text.length()-cstInputFileSuffix.length());
  text += ".nds";

  rein.open(text.c_str());
  message("trying to read '%s' ...\n", text.c_str());

  if (rein.fail())
  {
    throw cException("not able to open "+text, __FILE__, __LINE__);
    ExitApp();
  }
  else
  {
    trace("  file opened successfully ...");
  }

  // --------------------------------------------------------------------------
  //   read everything
  // --------------------------------------------------------------------------
  while (rein >> text)
  {
    // ------------------------------------------------------------------------
    //   look which room we're processing (1 or 2)
    // ------------------------------------------------------------------------
    if (text == "Seite")
    {
      rein >> seite;
      if ( (seite != 1) && (seite != 2) )
      {
        message("invalid value for 'seite' : %d\n", seite);
        throw cException("cMesh::parseNds() : invalid value for 'Seite'", __FILE__, __LINE__);
      }
    }
    // ------------------------------------------------------------------------
    //   ... nodes we have to evaluate ?
    // ------------------------------------------------------------------------
    else
    {
      // ------------------------------------------------------------------------
      //   looking for ':'
      // ------------------------------------------------------------------------
      int pos1 = (int)text.find(":");
      int pos2 = (int)text.rfind(":");

      // ------------------------------------------------------------------------
      //   ':' found ?
      // ------------------------------------------------------------------------
      if (pos1 >=0)
      {
        // ----------------------------------------------------------------------
        //   there is one ':' (pos1 == pos2) or are there two ':'
        // ----------------------------------------------------------------------
        if (pos1 == pos2)
        {
          std::string oben(text);
          std::string unten(text);

          oben.erase(0, pos1+1);
          unten.erase(pos1);


          for (int i=std::atoi(unten.c_str()); i<=std::atoi(oben.c_str()); i++)
            MyMesh.getNode( i )->setRoom( seite );
        }
        // ----------------------------------------------------------------------
        //   more than one ':' (pos1 != pos2)
        // ----------------------------------------------------------------------
        else
        {
          std::string oben(text);
          std::string unten(text);
          std::string step(text);

          oben.erase( 0, pos1+1 );
          oben.erase( oben.find(":") );
          unten.erase( pos1 );
          step.erase( 0, pos2+1 );

          for (int i=std::atoi(unten.c_str()); i<=std::atoi(oben.c_str()); i+=std::atoi(step.c_str()))
            MyMesh.getNode( i )->setRoom( seite );
        }
      }
      // ------------------------------------------------------------------------
      //   no ':' found
      // ------------------------------------------------------------------------
      else 
      {
        MyMesh.getNode( std::atoi(text.c_str()) )->setRoom( seite );
      }
    }
  }

  // --- close nds file
  rein.close();

  /*
  // -------------------------------------------------------------------------
  //   tell the fluid-structure coupling elements to which room they belong
  //   This information will be used later on to compute the sound intensity
  //   on each room's wall
  // -------------------------------------------------------------------------
  MyMesh.associateInterfaceElementsToRoom();
  */

  trace("finished.");
}


void cOutputLEV::computePressureLevel(const PetscReal &CurrentFrequency, cMesh &MyMesh)
{
  PetscInt FirstRow = 0; // first row of local portion of solution vector
  PetscInt LastRow = 0;  // last row of local portion of solution vector

  trace("  computing sound pressure level ... ");

  ierr = VecGetOwnershipRange(m_OutSolution[0], &FirstRow, &LastRow); INFAMCHKERRQ(ierr);

  // --- if the stream LevelFile is closed, this function is called
  //     for the first time. So we have to adjust the local range of
  //     rows and open the file.
  if ((m_GlobalNumRowsPerRoom[0] == 0) && (m_GlobalNumRowsPerRoom[1] == 0))
  {
    parseNdsFile(MyMesh);
    checkForSides(MyMesh);

    if (PetscGlobalRank == 0)
    {
      m_LevelFile.open( m_LevelFileName.c_str() );
      if (m_LevelFile.fail())
      {
        throw cException("Aborting! Unable to open " + m_LevelFileName, __FILE__, __LINE__);
      }

      // --- writing header of file
      m_LevelFile << "#--------------------------------------------------------" << std::endl;
      m_LevelFile << "#    f [Hz]        L1[dB]        L2[dB]      L1-L2[dB]   " << std::endl;
      m_LevelFile << "#                  Raum 1        Raum 2                  " << std::endl;
      m_LevelFile << "#--------------------------------------------------------" << std::endl;

      if (m_GlobalNumRowsPerRoom[0] <= 0)
        m_LevelFile << "#  no nodes specified for Room A !!" << std::endl;

      if (m_GlobalNumRowsPerRoom[1]  <= 0)
        m_LevelFile << "#  no nodes specified for Room B !!" << std::endl;


      m_LevelFile.setf(std::ios::scientific);
    }
  }


  // --- if the file is opened successfully, we'll now compute
  //     the sound pressure levels in each room ...
  trace("    computing local sum ... ");
  PetscScalar     value = 0.;

  PetscReal pEff;              // current effective value of the sound pressure 
  PetscReal pEffSumRoomA = 0.; // sum of effective values in room A
  PetscReal pEffSumRoomB = 0.; // sum of effective values in room B

  const PetscReal pRef           = 20.0e-6; // reference pressure in air
  const PetscReal numNodesRoomA  = (PetscReal)m_GlobalNumRowsPerRoom[0];
  const PetscReal numNodesRoomB  = (PetscReal)m_GlobalNumRowsPerRoom[1];


  short room; PetscInt row;
  for (ItMapNodes itN = MyMesh.getFirstNode(); itN != MyMesh.getLastNode(); itN++)
  {
    row = itN->second->getGlobalRow( fluid );

    if ((FirstRow <= row) && (row < LastRow))
    {
      room = itN->second->getRoom();
      VecGetValues(m_OutSolution[0], 1, &row, &value);

      //pEff = std::abs( value / M_SQRT2 );
#ifdef PETSC_USE_COMPLEX
      pEff = std::sqrt(value.real() * value.real() +  value.imag() * value.imag()) / M_SQRT2;
#else
      pEff = std::abs( value / M_SQRT2 );
#endif

      if (room == 1)
        pEffSumRoomA += pEff*pEff / (numNodesRoomA * pRef*pRef );
      else if (room == 2)
        pEffSumRoomB += pEff*pEff / (numNodesRoomB * pRef*pRef );
    }
  }


  // --- now we compute the sum over all processes and gather the result on rank 0
  trace("    computing global sum ... ");
  PetscReal globalSumRoomA = 0.;
  PetscReal globalSumRoomB = 0.;
  MPI_Reduce(&pEffSumRoomA, &globalSumRoomA, 1, MPIU_REAL, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&pEffSumRoomB, &globalSumRoomB, 1, MPIU_REAL, MPI_SUM, 0, PETSC_COMM_WORLD);


  PetscBool LaFlag  = PETSC_FALSE;
  PetscReal  GivenL1 = 0.;
  //PetscOptionsReal("-givenl1", "given SPL in sending room", "no manpage", 0.0, &GivenL1, &LaFlag);

  // --- output of the computed sound pressure level and the difference
  //     between them
  if (PetscGlobalRank == 0)
  {
    PetscReal L_A  = 0.;      // sound pressure level room A
    PetscReal L_B  = 0.;      // sound pressure level room B

    if (LaFlag == PETSC_TRUE) 
    {
      L_A = GivenL1;
    }
    else
    {
      if (globalSumRoomA > 0.)
        L_A = 10.* std::log10(globalSumRoomA);
      else
        trace("please check : LA <= 0.0");
    }

    if (globalSumRoomB > 0.)
      L_B = 10.* std::log10(globalSumRoomB);
    else
      trace("please check : LB <= 0.0");

    m_LevelFile << std::setw(15) << CurrentFrequency;
    m_LevelFile << std::setw(15) << L_A;
    m_LevelFile << std::setw(15) << L_B;
    m_LevelFile << std::setw(15) << L_A - L_B;
    m_LevelFile << std::endl << std::flush;
  }

  trace("  sound pressure level computed! ");
}



void cOutputLEV::writeMeshToFile(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    trace("no mesh output possible for LEV-Output ...");
  }
}

void cOutputLEV::writeMeshToFileGroup(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    trace("no mesh output for groups possible for LEV-Output ...");
  }
}

//! Creates a copy of current Object
cOutputfile* cOutputLEV::copy() {

  cOutputfile* copy_object = new cOutputLEV(*this);

  return copy_object;

}

std::ostream& cOutputLEV::write(std::ostream &os) const
{
  os << " want to compute pressure levels : " << wantOutput() << std::endl;

  return os;
}

std::istream& cOutputLEV::read(std::istream &is)
{
  return is;
}

std::ostream& cOutputLEV::writeXml(std::ostream &os) const
{
  os << "  <lev>" << wantOutput() << "</lev>" << std::endl;

  return os;
}
