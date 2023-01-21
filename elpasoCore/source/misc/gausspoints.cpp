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

#include "gausspoints.h"

cGaussPoints::cGaussPoints()
{
  // empty
}


cGaussPoints::~cGaussPoints()
{
  // empty
}


inline int cGaussPoints::Start(int level) const
{
  if (level <= 4)
  {
    double start = 0.0;
    start = ((double) level - 1.0) / 2.0 * (double) level;
    return (int) start;
  }
  else
  {
    std::cout<<" cGaussPoints::Start(int level) level must be < 4\n";
    ExitApp();
  }
  return -1;
}


cPoint cGaussPoints::getGaussPoint1D(int level, int nr) const
{
  cPoint gp;
  gp[0] = infam::gp4[Start(level) + nr];
  return gp;
}


PetscReal cGaussPoints::getGaussWeight1D(int level, int nr) const
{
  return (infam::gw4[Start(level)+nr]);
}


cPoint cGaussPoints::getGaussPoint3D(int level, int x, int y, int z) const
{
  cPoint gp;

  gp[0] = infam::gp4[Start(level) + x];
  gp[1] = infam::gp4[Start(level) + y];
  gp[2] = infam::gp4[Start(level) + z];

  return gp;
}

cPoint cGaussPoints::getGaussPoint3D(int level, int nr) const
{
  int r,x,y,z;

  r = nr % (level*level);
  z = (nr - r) / (level*level);
  r = (nr % (level*level)) % level;
  y = (nr % (level*level) - r) / level;
  x = (nr - y * level - z *level*level);

  return getGaussPoint3D(level, x, y, z);
}


cPoint cGaussPoints::getGaussPointTet3D(int level, int nr) const
{


  // -------------------------------------------------------------------------
  //   searching for starting position within vector
  // -------------------------------------------------------------------------


  int k, pos = 0;
  for (k=0;k<infam::anztetra;k++)
  {
    if (infam::anzgp[k] == level)
    {
      pos = infam::sttetra[k];
    }
  }

  pos += nr;
  pos *= 3;

  // -------------------------------------------------------------------------
  //   copy point's coordinates
  // -------------------------------------------------------------------------
  cPoint tp;
  tp[0] = infam::gp3t[pos];
  tp[1] = infam::gp3t[pos+1];
  tp[2] = infam::gp3t[pos+2];
  //std::cout<<"cGaussPoints::getGaussPointTet3D: "<<tp[0]<<" "<<tp[1]<<" "<<tp[2]<<"\n";
  return tp;

}

cPoint cGaussPoints::getGaussPointTet3DL(int level, int nr) const
{
  // -------------------------------------------------------------------------
  //   searching for starting position within vector
  // -------------------------------------------------------------------------
  int k, pos = 0;
  for (k=0;k<infam::anztetra;k++)
  {
    if (infam::anzgp[k] == level)
    {
      pos = infam::sttetra[k];
    }
  }
  //std::cout<<"Position1: "<<pos<<"\n";
  pos += 4*nr;
  //std::cout<<"Position2: "<<pos<<"\n";
  //pos *= 4;

  // -------------------------------------------------------------------------
  //   copy point's coordinates
  // -------------------------------------------------------------------------
  cPoint tp;
  tp[0] = infam::gp3t[pos];
  tp[1] = infam::gp3t[pos+1];
  tp[2] = infam::gp3t[pos+2];
  tp[3] = infam::gp3t[pos+3];
  //std::cout<<"cGaussPoints::getGaussPointTet3DL: "<<tp[0]<<" "<<tp[1]<<" "<<tp[2]<<" "<<tp[3]<<"\n";
  return tp;
}

cPoint cGaussPoints::getGaussPoint2D(int level, int x, int y) const
{
  cPoint gp;
  //std::cout<<"Start(level) + x="<<Start(level) + x<<"\n";
  //std::cout<<"Start(level) + y="<<Start(level) + y<<"\n";
  gp[0] = infam::gp4[Start(level) + x];
  gp[1] = infam::gp4[Start(level) + y];

  return gp;
}


cPoint cGaussPoints::getGaussPoint2D(int level, int nr) const
{
  //std::cout<<"cGaussPoints::getGaussPoint2D("<<level<<", int"<<nr<<") -->";
  int x,y;

  y = nr % level;
  x = (nr - y) / level;
  //std::cout<<"level="<<level<<", x="<<x<<", y="<<y<<"\n";
  return getGaussPoint2D(level, x, y);
}


PetscReal cGaussPoints::getGaussWeight3D(int level, int x, int y, int z) const
{

  return infam::gw4[Start(level) + x] * infam::gw4[Start(level) + y] * infam::gw4[Start(level) + z];
}


PetscReal cGaussPoints::getGaussWeight3D(int level, int nr) const
{
  int r,x,y,z;

  r = nr % (level*level);
  z = (nr - r) / (level*level);
  r = (nr % (level*level)) % level;
  y = (nr % (level*level) - r) / level;
  x = (nr - y * level - z *level*level);

  return getGaussWeight3D(level, x, y, z);
}


PetscReal cGaussPoints::getGaussWeightTet3D(int level, int nr) const
{

  // -------------------------------------------------------------------------
  //   searching for starting position within vector
  // -------------------------------------------------------------------------
  int k, pos = 0;
  for (k=0;k<infam::anztetra;k++)
  {
    if (infam::anzgp[k] == level)
    {
      pos = infam::sttetra[k];
    }
  }
  pos += nr;

  // -------------------------------------------------------------------------
  //   return weight
  // -------------------------------------------------------------------------
  //std::cout<<"cGaussPoints::getGaussWeightTet3D: "<<infam::gw3t[pos]<<"\n";
  return infam::gw3t[pos];


}

PetscReal cGaussPoints::getGaussWeightTet3DL(int level, int nr) const
{

  // -------------------------------------------------------------------------
  //   searching for starting position within vector
  // -------------------------------------------------------------------------
  int pos = 0;
  /*for (k=0;k<infam::anztetra;k++)
  {
  if (infam::anzgp[k] == level)
  {
  pos = infam::sttetra[k];
  }
  }*/
  switch(level)
  {
  case 6:
    pos=10;
    break;
  case 7:
    pos=11;
    break;
  case 8:
    pos=15;
    break;
  }

  pos += nr;

  // -------------------------------------------------------------------------
  //   return weight
  // -------------------------------------------------------------------------
  //std::cout<<"cGaussPoints::getGaussWeightTet3DL: "<<infam::gw3t[pos]<<"\n";
  return infam::gw3t[pos];

}

PetscReal cGaussPoints::getGaussWeight2D(int level, int x, int y) const
{
  return infam::gw4[Start(level) + x] * infam::gw4[Start(level) + y];
}


PetscReal cGaussPoints::getGaussWeight2D(int level, int nr) const
{
  int x,y;

  y = nr % level;
  x = (nr - y) / level;

  return getGaussWeight2D(level, x, y);
}


cPoint cGaussPoints::getGaussPointTria(int level, int nr) const
{
  // -------------------------------------------------------------------------
  //   searching for starting position within vector
  // -------------------------------------------------------------------------
  int k, pos = 0;
  for (k=0;k<infam::anztria;k++)
  {
    if (infam::arten[k] == level)
    {
      pos = infam::sttria[k];
    }
  }

  pos += nr;
  pos *= 3;

  // -------------------------------------------------------------------------
  //   copy point's coordinates
  // -------------------------------------------------------------------------
  cPoint p;
  p[0] = infam::gp3[pos];
  p[1] = infam::gp3[pos+1];
  p[2] = infam::gp3[pos+2];

  return p;
}

PetscReal cGaussPoints::getGaussWeightTria(int level, int nr) const
{
  // -------------------------------------------------------------------------
  //   searching for starting position within vector
  // -------------------------------------------------------------------------
  int k, pos = 0;
  for (k=0;k<infam::anztria;k++)
  {
    if (infam::arten[k] == level)
    {
      pos = infam::sttria[k];
    }
  }
  pos += nr;

  // -------------------------------------------------------------------------
  //   return weight
  // -------------------------------------------------------------------------
  return infam::gw3[pos];
}


