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

#include "elementfluid8.h"

cElementFluid8::cElementFluid8() : cElementFluid3d(8, 2) {
  // empty
}

cElementFluid8::cElementFluid8(const cElementFluid8 &other)
    : cElementFluid3d(other) {
  // empty
}

cElementFluid8::~cElementFluid8() {
  // empty
}

void cElementFluid8::getLocalCoordinatesOfElementsNodes(cMatrix &xyz) {
  // -- local coordinates of element's nodes
  xyz(0, 0) = -1.;
  xyz(0, 1) = -1.;
  xyz(0, 2) = -1.;  // node 1
  xyz(1, 0) = 1.;
  xyz(1, 1) = -1.;
  xyz(1, 2) = -1.;  // node 2 ...
  xyz(2, 0) = 1.;
  xyz(2, 1) = 1.;
  xyz(2, 2) = -1.;
  xyz(3, 0) = -1.;
  xyz(3, 1) = 1.;
  xyz(3, 2) = -1.;
  xyz(4, 0) = -1.;
  xyz(4, 1) = -1.;
  xyz(4, 2) = 1.;
  xyz(5, 0) = 1.;
  xyz(5, 1) = -1.;
  xyz(5, 2) = 1.;
  xyz(6, 0) = 1.;
  xyz(6, 1) = 1.;
  xyz(6, 2) = 1.;
  xyz(7, 0) = -1.;
  xyz(7, 1) = 1.;
  xyz(7, 2) = 1.;  // node 8
}

std::ostream &cElementFluid8::write(std::ostream &os) const {
  os << "Fluid8 element" << std::endl;
  os << "Id = " << getId() << std::endl;
  os << *m_Material << std::endl;

  for (int k = 0; k < getNumberOfNodes(); k++) os << *(m_Nodes[k]) << std::endl;

  return os;
}

std::ostream &cElementFluid8::writeXml(std::ostream &os) const {
  os << "<Fluid8>";
  os << "<Id>" << getId() << "</Id>";
  for (int k = 0; k < getNumberOfNodes(); k++)
    os << "<N>" << m_Nodes[k]->getId() << "</N>";
  os << "</Fluid8>";
  return os;
}

//! test by Meike
void cElementFluid8::computeVi(PetscReal &nimp, PetscScalar &Zn,
                               const Vec &FullSolution, const Vec &Stresses) {
  // if (m_ElementLoads.size() == 0)
  //   return;

#ifdef PETSC_USE_COMPLEX
  const int nnod = getNumberOfNodes();
  PetscInt index, ierr;
  PetscInt *indices = new PetscInt[3 * nnod];
  cElementVector gradP(3 * nnod);
  cElementVector valP(nnod);

  // --- get averaged grad(p) at each node
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalSeqId();

    indices[3 * k] = cstNumberOfStressDofs * index + gradpx;
    indices[3 * k + 1] = cstNumberOfStressDofs * index + gradpy;
    indices[3 * k + 2] = cstNumberOfStressDofs * index + gradpz;
  }

  ierr = VecGetValues(Stresses, 3 * nnod, indices, &(gradP[0]));
  INFAMCHKERRQ(ierr);
  delete[] indices;

  // --- extract fluid pressure from solution vector
  for (int k = 0; k < nnod; k++) {
    index = m_Nodes[k]->getGlobalRow(fluid);
    VecGetValues(FullSolution, 1, &index, &(valP[k]));
  }

  const PetscReal omega = m_Material->getOmega();
  const PetscReal rho = m_Material->getRho();

  std::complex<PetscReal> I(0., 1.);

  for (int k = 0; k < nnod; k++) {
    if (m_Nodes[k]->getId() == 115) {
      index = m_Nodes[k]->getGlobalSeqId();
      // particel velocity in z - Direktion !!only test case!!
      PetscScalar vn_ = -1. / (I * omega * rho) * gradP[3 * k + 2];
      // std::cout << vn_.real() << " " << vn_.imag() << std::endl;

      // std::cout << "index: " << index << " vi: "<< std::endl;
      /*std::cout.setf(std::ios::scientific);
      std::cout << std::setw(8) << m_Nodes[k]->getId();
      std::cout << std::setw(15) << omega;
      std::cout << std::setw(15) << Zn.real();
      std::cout << std::setw(15) << Zn.imag();

      //std::cout << " Pi: "<< std::endl;
      std::cout << std::setw(15) << valP[k].real();
      std::cout << std::setw(15) << valP[k].imag();*/

      // std::cout << " Zi: "<< std::endl;
      PetscScalar Zi = valP[k] / vn_;
      /*std::cout << std::setw(15) << Zi.real();
      std::cout << std::setw(15) << Zi.imag() << std::endl;*/
      // std::cout << "Point: " << *(m_Nodes[k]) << " Z meike real " <<
      // Zi.real() <<  " imag " << Zi.imag() << std::endl;
      std::cout << Zi.real() << " " << Zi.imag() << std::endl;

      PetscScalar Z0 = 1.2 * 346.0;
      PetscReal absor = 1 - PetscAbsScalar((Zi - Z0) / (Zi + Z0)) *
                                PetscAbsScalar((Zi - Z0) / (Zi + Z0));
      // std::cout << absor << std::endl;
    }
  }

#else
  trace("********************************************************");
  trace("* computeVi is only available in frequency domain      *");
  trace("* Element isn't suitable for time domain computations. *");
  trace("* Please rebuild application using complex numbers     *");
  trace("********************************************************");
#endif
}
