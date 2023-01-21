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

#ifndef INFAM_BASICSDOF_H
#define INFAM_BASICSDOF_H

//! degrees of freedom used in the program
const short cstNumberOfKnownDofs = 21;

//! @brief synonyms for the degrees of freedom
//! they have to be numbered sequentially because they
//! may be used within loops
enum eKnownDofs {
  disp_x1 = 0, /*!< displacement in x-direction */
  disp_x2 = 1, /*!< displacement in y-direction */
  disp_x3 = 2, /*!< displacement in z-direction */
  disp_w1 = 3, /*!< rotation around x-axis */
  disp_w2 = 4, /*!< rotation around y-axis */
  disp_w3 = 5, /*!< rotation around z-axis */
  pore0 = 6,   /*!< pore pressure (constant term) */
  pore1 = 7,   /*!< pore pressure (linear term) */
  pore2 = 8,   /*!< pore pressure (quadratic term) */
  pore3 = 9,   /*!< pore pressure (cubic term) */
  disp_xd3 = 10,
  disp_wd1 = 11,
  disp_wd2 = 12,
  disp_dwdx = 13,  /*!< w,x (Bernoulli beam/Kirchhoff plate */
  disp_dwdy = 14,  /*!< w,y (Bernoulli beam/Kirchhoff plate */
  disp_dwdxy = 15, /*!< w,xy (Bernoulli beam/Kirchhoff plate */
  fluid = 16,      /*!< fluid pressure */
  disp_z_1 =
      17, /*!< displacement in z-direction (linear term) - (PoroDiscKienzler)*/
  disp_z_3 =
      18, /*!< displacement in z-direction (cubic term) - (PoroDiscKienzler)*/
  disp_x1_2 = 19, /*!< displacement in x-direction (quadratic term) -
                     (PoroDiscKienzler)*/
  disp_x2_2 = 20  /*!< displacement in y-direction (quadratic term) -
                     (PoroDiscKienzler)*/
};

const eKnownDofs cstAllDofs[] = {
    disp_x1,   disp_x2,    disp_x3, disp_w1,  disp_w2,  disp_w3,   pore0,
    pore1,     pore2,      pore3,   disp_xd3, disp_wd1, disp_wd2,  disp_dwdx,
    disp_dwdy, disp_dwdxy, fluid,   disp_z_1, disp_z_3, disp_x1_2, disp_x2_2};

//! the degrees of freedom as text - used for some output
const std::string sKnownDofsIdentifiers[] = {
    "disp_x1",  "disp_x2",   "disp_x3",   "disp_w1",    "disp_w2",  "disp_w3",
    "pore0",    "pore1",     "pore2",     "pore3",      "disp_xd3", "disp_wd1",
    "disp_wd2", "disp_dwdx", "disp_dwdy", "disp_dwdxy", "fluid",    "disp_z_1",
    "disp_z_3", "disp_x1_2", "disp_x2_2"};

#endif
