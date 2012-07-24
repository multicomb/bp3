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

#ifndef BP3MATRIX_H 
#define BP3MATRIX_H 1

#include "misc.h"

class Matrix
{
  // This is for real-symmetric-square matrices only.
  // However, I have tried to make a transition to a general matrix easy.
  // For reasons of efficiency, I don't do any checking of array bounds...
  
private:
  vector<double> data;
  unsigned row, col, si;
  
public:
  Matrix(const unsigned i = 2)
  { si = i*i; row = col = i; data.resize(si);                                     }
  Matrix(const unsigned i, const double *init)
  { si = i*i; row = col = i; data.resize(si);
    for (unsigned i = 0; i < si; i++) data[i] = init[i];                          }
  ~Matrix()                                        { ;                            }
  void resize(unsigned i)
  { data.resize(i*i); row = col = i; si = i*i;                                    }
  void resize(unsigned i, unsigned j, bool eigen = true)
  { data.resize(i*j); row = i; col = j; si = i*j;                                 }
  unsigned size()  const                           { return si;                   }
  unsigned rsize() const                           { return row;                  }      
  unsigned csize() const                           { return col;                  }
  double *array()                                  { return &data[0];             }
  double vec(const unsigned i) const               { return data[i];              }
  double  operator()(unsigned i, unsigned j) const { return data[i+col*j];        }
  double &operator()(unsigned i, unsigned j)       { return data[i+col*j];        }
  Matrix &assign(const double *in)
  {for (unsigned i = 0; i < si; i++) data[i] = in[i]; return *this; }
  Matrix &assign(const double *in, const unsigned len)
  {for (unsigned i = 0; i < len; i++) data[i] = in[i]; resize(len); return *this; }
  /*
  Matrix& operator*=(const double c)
  {for (unsigned i = 0; i < si; i++) data[i]     *= c; return *this;         }
  Matrix& operator/=(const double c)
  {for (unsigned i = 0; i < si; i++) data[i]     /= c; return *this;         }
  Matrix& operator*(const Matrix &a, const Matrix &b)
  {
    Matrix c(a.rsize(); b.csize());
    for (unsigned i       = 0; i < a.rsize(); i++)
      for (unsigned j     = 0; j < b.csize(); j++)
	c[i+b.csize()*j] += 

    return c;
  }
  */
};

#endif
