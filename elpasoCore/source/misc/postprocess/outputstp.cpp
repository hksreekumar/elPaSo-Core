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

#include "outputstp.h"

cOutputSTP::cOutputSTP()
{  
}

cOutputSTP::cOutputSTP(cOutputSTP &other):
cOutputfile(other)
{
}

cOutputSTP::~cOutputSTP()
{
}


void cOutputSTP::writeResultsToFile(cMesh &MyMesh, const int &NrStep, const double &Step)
{
  if (wantOutput() > 0 || checkForFilter() == true) {
    if(filter_init == false)
      initializeFilter(MyMesh);

    writeResultSTP(NrStep, Step, MyMesh.getFilename());
  }
}

void cOutputSTP::writeResultSTP(const int &step, const double &frequenz, const std::string &Filename)
{
  // --- check if the user really wants to have output
  if ((wantOutput() == 0) && (checkForFilter() != true)) {
    return;
  }

  // --- collect solution on rank 0
  Vec        FullVector;
  PetscInt   GlobalSize;
  VecScatter ctx;

  VecScatterCreateToZero(m_OutSolution[0], &ctx, &FullVector);
  VecScatterBegin(ctx,m_OutSolution[0],FullVector,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,m_OutSolution[0],FullVector,INSERT_VALUES,SCATTER_FORWARD);

  VecGetSize(FullVector, &GlobalSize);


  // --- filter some nodal results from solution vector
  if (checkForFilter() == true) {
    filterFromSolution(frequenz, FullVector);
  }
  if (wantOutput() == 0) {
    return;
  }

  // --- full output of solution vector
  if (wantOutput() > 0) {
    ogzstream         raus;
    std::stringstream text;
    std::string       CurrentName(Filename);

    // --- compose current filename
    if (CurrentName.rfind(cstInputFileSuffix) == (CurrentName.length()-cstInputFileSuffix.length()))
      CurrentName.erase(CurrentName.length()-cstInputFileSuffix.length());

    text << ".";
    text.width(8);
    text.fill('0');
    text << step << ".stp.gz";
    CurrentName += text.str();
    text.str("");

    message("  writing %s ...\n", CurrentName.c_str());

    if (PetscGlobalRank == 0)
    {
      PetscScalar SingleValue;

      raus.open(CurrentName.c_str());
      if (raus.fail())
      {
        throw cException("Aborting! not able to create " + CurrentName, __FILE__, __LINE__);
      }

      raus.setf(std::ios::scientific);
      raus << frequenz << std::endl;

      // --- write portion of solution located on rank 0
      for (PetscInt k=0; k<GlobalSize; k++)
      {
        VecGetValues(FullVector, 1, &k, &SingleValue);
#ifdef PETSC_USE_COMPLEX
        raus << std::setw(15) << SingleValue.real();
        raus << std::setw(15) << SingleValue.imag();
#else
        raus << std::setw(15) << SingleValue;
#endif
        raus << std::endl;
      }
    }

    if (PetscGlobalRank == 0)
      raus.close();
  }


  VecScatterDestroy(&ctx);
  VecDestroy(&FullVector);

  trace("  output finished.");
}


void cOutputSTP::writeMeshToFile(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    trace("no mesh output possible for STP-Output ...");
  }
}

void cOutputSTP::writeMeshToFileGroup(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    trace("no mesh output for groups possible for STP-Output ...");
  }
}

//! Creates a copy of current Object
cOutputfile* cOutputSTP::copy() {

  cOutputfile* copy_object = new cOutputSTP(*this);

  return copy_object;

}

std::ostream& cOutputSTP::write(std::ostream &os) const
{
  os << " want to write .stp.gz files     : " << wantOutput() << std::endl;

  return os;
}


std::istream& cOutputSTP::read(std::istream &is)
{
  return is;
}


std::ostream& cOutputSTP::writeXml(std::ostream &os) const
{
  os << "  <stp>" << wantOutput() << "</stp>" << std::endl;

  return os;
}
