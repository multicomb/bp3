/*
  Copyright (c) Navraj Pannu       1998 - 2005
  Copyright (c) Leiden University  2002 - 2005

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.     
*/

#ifndef BP3INPUT_H
#define BP3INPUT_H	1

#include "bp3likelihood.h"
#include "crystal.h"
#include "input.h"
#include "model.h"
#include "minimizer.h"

class Bp3input : public Input
{
 public:
  Bp3input(Model &, Crystal &, Bp3likelihood &, Minimizer &);
  ~Bp3input(){;}
  
 private:
  Bp3likelihood &like;
  
  void CCP4parse();
  bool bp3keywords();
  void refikeyword() const;
  void phaskeyword() const;
  void checkeyword() const;
  void handkeyword() const;
  void diffkeyword() const;
  void shelkeyword() const;
  void fhoukeyword() const;
  void threkeyword() const;
  void check();
  void print() const
    { printf("\nAtomic parameters inputted\n"); mdl.print(true);}  
};

#endif
