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

#ifndef MISC_H
#define MISC_H	1

#include <functional>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "cmtzlib.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

// working precision of all floating point variables
// typedef double work;
typedef float work;

// Numerical Constants                 
const double SMALLESTD     = 1.0e-100;
const double DSMALL        = 1.0e-15;
const work   WEPSILON      = 1.0e-5;
const double EPSILON       = 1.0e-5;
const double ZERO          = 0.0;
const double QUARTER       = 0.25;
const double THIRD         = 0.33333333333333333333333;
const double HALF          = 0.5;
const double TWOTHIRD      = 0.66666666666666666666666;
const double ONE           = 1.0;
const double TWO           = 2.0;
const double THREE         = 3.0;
const double FOUR          = 4.0;
const double FIVE          = 5.0;
const double TEN           = 10.0;
const double TWENTYFIVE    = 25.0;
const double NINETY        = 90.0;
const double HUNDREDTWENTY = 120.0;
const double THOUSAND      = 1000.0;
const double PI            = 3.14159265358979323846264;
const double TWOPI         = 6.28318530717958623199593;
const double TWOPIINV      = ONE/6.28318530717958623199593;
const double TWOPI2        = 19.7392088021787172376690;
const double EIGHTPI2      = 78.9568352087148689506759;
const double DEGREEtoRAD   = 0.01745329251994329576924;

// number of carbons/oxygens etc per residue/nucleotide
const double CARBONRES     = 4.869;
const double NITRORES      = 1.351;
const double OXYGENRES     = 1.351;
const double SULFURRES     = 0.051;
const double HYDRORES      = 8.000;
const double CARBONNUC     = 9.683;
const double NITRONUC      = 1.351;
const double OXYGENNUC     = 6.424;
const double PHOSNUC       = 1.000;
const double HYDRONUC      = 12.000;

// minimum eigenvalue allowed for positive definite matrices
const double MINEIG        = 0.00001;

// bp3's "not used" flag for diffraction data
const work NOTUSED         = -2525252525.25;

void CCP4Banner(const int argc, char **argv);
void Termination();
void Bp3Error(const string, const string, const bool = true); 
void Bp3Warning(const string, const string); 
void Bp3Result(const string); 

// mathematical functions
double i0e(double, int, int);                      // exp(-x) I0(x)
double i1e(double, int, int);                      // exp(-x) I1(x)
double chbevl(double, void *, int, int);
double sim(const double x, int, int);             // I1(x)/I0(x)
double simoverx(const double x, int, int);        // I1(x)/(I0(x)*x)

// for window's compatibility, atanh has to be defined
// thanks Martyn/CCP4!

#if defined _MVS && !defined atanh
#define atanh(x) 0.5*log((1+x)/(1-x));
#endif

template< class T >
const T &between(const T &a, const T &b, const T &c)
{
  return std::max(a, std::min(b, c));
}

template< class T >
const T interpolate(const T &x, const T &a, const T &b, const T &ya, const T &yb)
{
  // linear interpolation

  const T c   = b - a;
  if (fabs(c) > ZERO)
    return ((b - x)*ya + (x - a)*yb)/c;
  else
    return yb;
}

template< class T >
const T dsign(const T &a, const T &b)
{
  return  ( (b >= 0) ? fabs(a) : -fabs(a) );
}

template< class T >
const T sign(const T &a)
{
  return  ( (a >= 0) ? (T) ONE : -(T) ONE );
}

template< class T >
const T dot_prod(const vector<T> &a, const vector<T> &b)
{
  if (a.size()   != b.size())
    Bp3Error("dot_prod","vectors are a different length");

  T dot(0);

  for (unsigned i = 0; i < a.size(); i++)
    dot          += a[i]*b[i];
  
  return dot;
}

template< class T >
const T diff(const vector<T> &a, const vector<T> &b)
{
  if (a.size()   != b.size())
    Bp3Error("diff","vectors are a different length");

  T diff(0);

  for (unsigned i = 0; i < a.size(); i++)
    diff         += fabs(a[i] - b[i]);
  
  return diff;
}

template< class T >
const T norm(const vector<T> &a)
{
  T norm(ZERO);

  for (unsigned i = 0; i < a.size(); i++)
    norm         += a[i]*a[i];
  
  return sqrt(norm);
}

#endif
