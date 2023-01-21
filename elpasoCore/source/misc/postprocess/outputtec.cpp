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

#include "outputtec.h"

cOutputTEC::cOutputTEC()
{  
}

cOutputTEC::cOutputTEC(cOutputTEC &other):
cOutputfile(other)
{
}


cOutputTEC::~cOutputTEC()
{
}


void cOutputTEC::writeResultsToFile(cMesh &MyMesh, const int &NrStep, const double &Step)
{
  if (wantOutput() > 0) {

    exportSolutionToTecplot(&MyMesh, NrStep, Step);
  }
}




void cOutputTEC::exportSingleTypeToTecPlot(const eElementShape &Shape, cMesh &MyMesh, 
                                           std::ostream &output, const PetscReal &step, ePhysicsType PhysicsType,
                                           const int &GroupId)
{
  std::string  TecPlotElementType; // TecPlot's identifier for an elementtype
  std::string  TecPlotZoneName;    // name for zone that will appear in TecPlot's frame
  int          NodesNeeded = 0;    // number of nodes of specific type that build the element

  std::set<PetscInt> NodeIds;      // ids of the nodes for elements of specific type
  std::set<PetscInt> ElementIds;   // ids of elements for a specific type

  switch(Shape)
  {
  case Beam:
    TecPlotElementType = "LINESEG";
    TecPlotZoneName    = "Beams";
    NodesNeeded        = 2;
    break;
  case Tria:
    TecPlotElementType = "TRIANGLE";
    TecPlotZoneName    = "Triangles";
    NodesNeeded        = 3;
    break;
  case Quadrilateral:
    TecPlotElementType = "QUADRILATERAL";
    TecPlotZoneName    = "Quadrilaterals";
    NodesNeeded        = 4;
    break;
  case Tetrahedron:
    TecPlotElementType = "TETRAHEDRON";
    TecPlotZoneName    = "Tetrahedrons";
    NodesNeeded        = 4;
    break;
  case Hexahedron:
    TecPlotElementType = "BRICK";
    TecPlotZoneName    = "Hexahedrons";
    NodesNeeded        = 8;
    break;
  case QuadSerendipity:
    trace(" warning: serendipity elements will not be");
    trace("          exported to Tecplots DAT file !!");
    return;
    break;
  }

  message("  Output of '%s'\n", TecPlotZoneName.c_str());

  if (GroupId == -1) {
    // --- loop over elements and look for the specified shape
    //     If this shape is found read as much nodes as needed to draw the element.
    //     The mid-side-nodes of quadratic elements are ignored.
    for (ItMapElements it = MyMesh.getFirstElement(); it != MyMesh.getLastElement(); it++) {
      if (it->second->getElementShape() == Shape) {
        // first is used for export only the mesh, second one if we want to export
        // the solution at the nodes
        if ((PhysicsType == 0) || (it->second->getPhysicsType() == PhysicsType)) {
          ElementIds.insert(it->second->getId());

          for (int k=0; k<NodesNeeded; k++) {
            NodeIds.insert( it->second->getNode(k)->getId() );
          }
        }
      }
    }
  }
  else {
    cElementFEM *ptrElement;
    std::set<PetscInt>::iterator itSet;

    PetscPrintf(PETSC_COMM_SELF, "Gruppe : %d \n", GroupId);

    for (ItMapGroups itG = MyMesh.getFirstGroup(); itG != MyMesh.getLastGroup(); itG++) 
    {
      if (itG->second->getId() == GroupId) 
      {
        for (itSet = itG->second->getFirstElementId(); itSet != itG->second->getLastElementId(); itSet++) 
        {
          ptrElement = MyMesh.getElement( *itSet );

          if (((ptrElement->getElementShape() == Shape) && (PhysicsType == 0)) ||
            (ptrElement->getPhysicsType() == PhysicsType)) 
          {
            ElementIds.insert( *itSet );
            for (int k=0; k<ptrElement->getNumberOfNodes(); k++) 
            {
              NodeIds.insert( ptrElement->getNode(k)->getId() );
            }
          }

          ptrElement = NULL;
        }
      }
    }
  }


  // --- leave if no appropriate elements are found
  //     (you need to check on all ranks)
  PetscInt GlobalSum = 0;
  PetscInt LocalSum = ElementIds.size();

  MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);
  if (GlobalSum == 0) {
    return;
  }


  // --- collect all nodes needed to display the
  //     mesh on rank 0. First, the number of nodes found
  //     on all procs is collected on rank 0
  PetscInt LocalSize = NodeIds.size();
  std::vector<PetscInt> NodeIdsPerProc(PetscGlobalSize, 0);
  MPI_Gather( &LocalSize, 1, MPIU_INT, &(NodeIdsPerProc[0]), 1, MPIU_INT, 0, PETSC_COMM_WORLD);

  // --- now all nodes are sent to rank 0 and added to the set of nodes
  for (int rank=1; rank<PetscGlobalSize; rank++)
  {
    if (PetscGlobalRank == rank)
    {
      std::vector<PetscInt> buffer(NodeIds.size());
      std::copy(NodeIds.begin(), NodeIds.end(), buffer.begin());
      MPI_Send( &(buffer[0]), buffer.size(), MPIU_INT, 0, rank, PETSC_COMM_WORLD );
    }
    else if (PetscGlobalRank == 0)
    {
      MPI_Status status;
      std::vector<PetscInt> buffer( NodeIdsPerProc[rank], 0 );
      MPI_Recv( &(buffer[0]), buffer.size(), MPIU_INT, rank, rank, PETSC_COMM_WORLD, &status );

      for (int k=0; k<NodeIdsPerProc[rank]; k++)
      {
        NodeIds.insert( buffer[k] );
      }
    }
  }


  // --- collect solution on rank 0
  Vec        fullSolution;
  VecScatter ctx;

  if(m_OutSolution!=NULL)
  {
    VecScatterCreateToZero(m_OutSolution[0], &ctx, &fullSolution);
    VecScatterBegin(ctx,m_OutSolution[0],fullSolution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,m_OutSolution[0],fullSolution,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx);
  }

  // --- writing header of current zone
  std::set<PetscInt>::const_iterator itSet;

  if (PetscGlobalRank == 0) {
    std::stringstream Temp;
    Temp << step;
    output << "TITLE = \"" << TecPlotZoneName << " " << Temp.str() << "\" " << std::endl;

    std::string Variables = "VARIABLES = \"X\" \"Y\" \"Z\"";

    if (PhysicsType == Acoustics)
      Variables += " \"magnitude\" \"Log10p\" \"Lp dB\" \"anglephi\" ";
    else if (PhysicsType == Elastodynamics)
      Variables += " \"disp_x\" \"disp_y\" \"disp_z\" ";

    output << Variables << std::endl;

    output << "ZONE ";
    output << "N=" << NodeIds.size() << ", ";
    output << "E=" << GlobalSum << ", ";
    output << "F=FEPOINT, ET=" << TecPlotElementType << std::endl;
    output.setf(std::ios::scientific);

    // --- output of nodal coordinates
    //     we are only printing the nodes we've found
    //     in the sought for elements
    cNode *ptrNode = NULL;

    for (itSet=NodeIds.begin(); itSet!=NodeIds.end(); itSet++)
    {
      ptrNode = MyMesh.getNode( *itSet );

      // --- output of node's coordinates
      for (int k=0; k<3; k++)
      {
        PetscReal xi = (*ptrNode)[k];

        // --- may consider computed displacements
        if (PhysicsType == Elastodynamics)
        {
          PetscInt    row;
          PetscScalar Value = 0.;

          if (ptrNode->checkIfActive( k ) == true)
          {
            row = ptrNode->getGlobalRow( k );
            VecGetValues(fullSolution, 1, &row, &Value);

#ifdef PETSC_USE_COMPLEX
            xi += Value.real();
#else
            xi += Value;
#endif
          }
        }
        output << std::showpoint << xi << " ";
      }

      if (m_OutSolution != NULL)
      {
        if (PhysicsType == Acoustics)
        {
          PetscInt    row = ptrNode->getGlobalRow(fluid);
          PetscScalar Value;
          VecGetValues(fullSolution, 1, &row, &Value);
          PetscReal Magnitude = std::abs(Value);
          PetscReal LpdB      = 20. * std::log10(Magnitude*50000.);
          PetscReal logLp     = std::log10(Magnitude);

          output << Magnitude << " " << logLp << " " << LpdB << " ";
#ifdef PETSC_USE_COMPLEX
          output << computePhi(Value);
#else
          output << 0.;
#endif
        }
        else if(PhysicsType == Elastodynamics)
        {
          PetscScalar delta_x=0.0;
          PetscScalar delta_y=0.0;
          PetscScalar delta_z=0.0;
          PetscInt  row;

          //Output of abs values in case of complex calculation
#ifdef PETSC_USE_COMPLEX
          PetscReal dx=0.0;
          PetscReal dy=0.0;
          PetscReal dz=0.0;
#endif
          if(ptrNode->checkIfActive( disp_x1 ) == true)
          {
            row = ptrNode->getGlobalRow(disp_x1);
            VecGetValues(fullSolution, 1, &row, &delta_x);
          }
#ifdef PETSC_USE_COMPLEX
          //Output of abs values in case of complex calculation
          //delta_x=0.0;
          dx=abs(delta_x);
#endif
          if(ptrNode->checkIfActive( disp_x2 ) == true)
          {
            row = ptrNode->getGlobalRow(disp_x2);
            VecGetValues(fullSolution, 1, &row, &delta_y);
          }
#ifdef PETSC_USE_COMPLEX
          //delta_y=0.0;
          //Output of abs values in case of complex calculation
          dy=abs(delta_y);
#endif		
          if(ptrNode->checkIfActive( disp_x3 ) == true)
          {
            row = ptrNode->getGlobalRow(disp_x3);
            VecGetValues(fullSolution, 1, &row, &delta_z);
          }
#ifdef PETSC_USE_COMPLEX
          //delta_z=0.0;
          //Output of abs values in case of complex calculation
          dz=abs(delta_z);
#endif		

#ifdef PETSC_USE_COMPLEX
          //Output of abs values in case of complex calculation
          output << dx << " " << dy << " " << dz ;
#else
          output << delta_x << " " << delta_y << " " << delta_z ;
#endif		
        }
        else {
          // --- the number of zeros corresponds to the
          //     number of values computed for the acoustic domain
          output << "0. 0. 0. 0.";
        }
      }

      output << std::endl;
    }
    output.unsetf(std::ios::scientific);

    ptrNode = NULL;
  } // end rank 0 only

  if(m_OutSolution!=NULL)
  {
    VecDestroy(&fullSolution);
  }

  // --- copy node ids into a vector, first only on
  //     rank 0. In order to get correct connectivity
  //     tables the vector must be sent to all other
  //     ranks.
  PetscInt VecSizeRank0 = NodeIds.size(); // correct value currently only on rank 0 !!
  MPI_Bcast ( &VecSizeRank0, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
  std::vector<PetscInt> VecNodeIds( VecSizeRank0 );

  if (PetscGlobalRank == 0)
    std::copy(NodeIds.begin(), NodeIds.end(), VecNodeIds.begin());

  MPI_Bcast ( &(VecNodeIds[0]), VecSizeRank0, MPIU_INT, 0, PETSC_COMM_WORLD);


  // --- output of the elements
  cElement *ptrElement = NULL;
  PetscInt  puffer[8];

  if (PetscGlobalRank == 0)
  {
    for (itSet=ElementIds.begin(); itSet!=ElementIds.end(); itSet++)
    {
      ptrElement = MyMesh.getElement( *itSet );
      for (int k=0; k<NodesNeeded; k++)
      {
        // --- now compute the index of the value within the vector
        //     thx to de.comp.lang.iso-c++ :-)
        output << (int)(std::find(VecNodeIds.begin(), VecNodeIds.end(), ptrElement->getNode(k)->getId()) - VecNodeIds.begin()) + 1 << " ";
      }
      output << std::endl;
    }
    ptrElement = NULL;
  }
  else {
    // --- all other ranks send their elements
    int count = 0;
    for (itSet=ElementIds.begin(); itSet!=ElementIds.end(); itSet++)
    {
      ptrElement = MyMesh.getElement( *itSet );

      // --- zero entries
      for (int k=0; k<8; k++)
        puffer[k] = 0;

      // --- determine node ids
      for (int k=0; k<NodesNeeded; k++)
      {
        puffer[k] = (int)(std::find(VecNodeIds.begin(), VecNodeIds.end(),
          ptrElement->getNode(k)->getId()) - VecNodeIds.begin()) + 1;
      }

      ptrElement = NULL;

      // --- send data
      //     blocking send -> second block will be sent if rank 0 fetches the first one
      MPI_Send( puffer, 8, MPIU_INT, 0, 1000+count, PETSC_COMM_WORLD );
      count++;
    }

    // --- send termination sign (first value is -1)
    puffer[0] = -1;
    MPI_Send( puffer, 8, MPIU_INT, 0, 1000+count, PETSC_COMM_WORLD );
  }


  // --- receive data from other ranks
  if (PetscGlobalRank == 0)
  {
    MPI_Status status;

    PetscInt rank = 1;
    PetscInt count = 0;

    while (rank<PetscGlobalSize)
    {
      MPI_Recv( puffer, 8, MPIU_INT, rank, 1000+count, PETSC_COMM_WORLD, &status );
      count++;
      if (puffer[0] == -1) {
        rank++;
        count = 0;
      }
      else {
        for (int k=0; k<NodesNeeded; k++)
          output << puffer[k] << " ";
        output << std::endl;
      }
    }
  }
}


void cOutputTEC::exportMeshToTecplot(cMesh *myMesh, bool byGroup)
{
  std::ofstream Output;
  std::string   FilenameRoot;
  std::string   Filename;

  // --- join filename of outputfile
  createFilenameRoot(FilenameRoot, myMesh);
  Filename = FilenameRoot + ".dat";

  // -- open Tecplot DAT file
  Output.open(Filename.c_str());


  if(byGroup == false)
  {
    trace("export mesh to Tecplot ... ");
    exportSingleTypeToTecPlot(Beam, *myMesh, Output);
    exportSingleTypeToTecPlot(Tria, *myMesh, Output);
    exportSingleTypeToTecPlot(Quadrilateral, *myMesh, Output);
    exportSingleTypeToTecPlot(Tetrahedron, *myMesh, Output);
    exportSingleTypeToTecPlot(Hexahedron, *myMesh, Output);
  }
  else
  {
    // -- if different element types exist for a single group multiple zones are
    //    exported for a single group ...
    trace("export mesh to Tecplot ordered by group ...");
    for (ItMapGroups itG = myMesh->getFirstGroup(); itG != myMesh->getLastGroup(); itG++) {
      exportSingleTypeToTecPlot(Beam, *myMesh, Output, 0.0, Undefined, itG->second->getId());
      exportSingleTypeToTecPlot(Tria, *myMesh, Output, 0.0, Undefined, itG->second->getId());
      exportSingleTypeToTecPlot(Quadrilateral, *myMesh, Output, 0.0, Undefined, itG->second->getId());
      exportSingleTypeToTecPlot(Tetrahedron, *myMesh, Output, 0.0, Undefined, itG->second->getId());
      exportSingleTypeToTecPlot(Hexahedron, *myMesh, Output, 0.0, Undefined, itG->second->getId());
    }
  }

  // --- close file
  Output.close();
}

void cOutputTEC::exportSolutionToTecplot(cMesh *myMesh, const PetscInt &NumStep, const PetscReal &Step)
{
  std::string   Filename;
  std::string   FilenameRoot;
  std::stringstream Temp;

  trace("export solution to TecPlot ...");

  if (PetscGlobalRank == 0)
  {
    // --- join filename of outputfile
    createFilenameRoot(FilenameRoot, myMesh);

    Temp << ".";
    Temp.width(8);
    Temp.fill('0');
    Temp << NumStep << ".dat";
    Filename += FilenameRoot + Temp.str();

    // -- open Tecplot DAT file
    TecplotOutput.open(Filename.c_str());
  }

  exportSingleTypeToTecPlot(Beam,          *myMesh, TecplotOutput, Step, Elastodynamics);
  exportSingleTypeToTecPlot(Tria,          *myMesh, TecplotOutput, Step, Elastodynamics);
  exportSingleTypeToTecPlot(Quadrilateral, *myMesh, TecplotOutput, Step, Elastodynamics);
  exportSingleTypeToTecPlot(Hexahedron   , *myMesh, TecplotOutput, Step, Elastodynamics);
  exportSingleTypeToTecPlot(Hexahedron   , *myMesh, TecplotOutput, Step, Acoustics);
  exportSingleTypeToTecPlot(Quadrilateral, *myMesh, TecplotOutput, Step, Acoustics);

  // --- close file
  if (PetscGlobalRank == 0)
    TecplotOutput.close();
}


PetscReal cOutputTEC::computePhi(std::complex<PetscReal> &number) const
{
  if (number.real() > 0.)
    return (std::atan(number.imag() / number.real()));
  else if ((number.real() == 0.) && (number.imag() > 0.))
    return (M_PI / 2.);
  else if ((number.real() == 0.) && (number.imag() < 0.))
    return (-1. * M_PI / 2.);
  else if ((number.real() < 0.) && (number.imag() >= 0.))
    return (atan(number.imag() / number.real()) + M_PI);
  else if ((number.real() < 0.) && (number.imag() < 0.))
    return (atan(number.imag() / number.real()) - M_PI);
  else
    return 0.;
}




void cOutputTEC::writeMeshToFile(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    exportMeshToTecplot(&MyMesh, false);
  }
}

void cOutputTEC::writeMeshToFileGroup(cMesh &MyMesh)
{
  if (wantOutput() > 0) {

    exportMeshToTecplot(&MyMesh, true);
  }
}

//! Creates a copy of current Object
cOutputfile* cOutputTEC::copy() {

  cOutputfile* copy_object = new cOutputTEC(*this);

  return copy_object;

}

std::ostream& cOutputTEC::write(std::ostream &os) const
{
  os << " want a TecPlot file             : " << wantOutput() << std::endl;

  return os;
}


std::istream& cOutputTEC::read(std::istream &is)
{
  return is;
}


std::ostream& cOutputTEC::writeXml(std::ostream &os) const
{
  os << "  <tec>" << wantOutput() << "</tec>" << std::endl;

  return os;
}
