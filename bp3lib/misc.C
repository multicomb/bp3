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

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include "misc.h"
#include "ccp4_general.h"
#include "ccp4_program.h"
#include "bp3_version.h"

void CCP4Banner(const int argc, char **argv)
{
  CCP4::ccp4fyp(argc,argv);
  CCP4::ccp4_prog_vers(version.c_str());
  CCP4::ccp4_banner();
}

void Termination()
{
  CCP4::ccperror(0,"Normal termination"); 
}

void Bp3Error(const string methodname, const string errormessage, const bool error)
{
  // Reports a message from a method where an                      
  // error has occurred and terminates the program                   

  printf("<B><FONT COLOR=\"#FF0000\"><!--SUMMARY_BEGIN-->\n");
  printf("\nBPIII: Error in method:  %s  **%s**\n\n",
	 methodname.c_str(), errormessage.c_str());
  printf("<!--SUMMARY_END--></FONT></B>\n");

  if (error)  
    exit(1);
}

void Bp3Warning(const string methodname, const string warningmessage)
{
  // Reports a warning message from a method
  printf("<B><FONT COLOR=\"#FF0000\"><!--SUMMARY_BEGIN-->\n");
  printf("\nBPIII: Warning in method:  %s  **%s**\n\n",
	 methodname.c_str(), warningmessage.c_str());
  printf("$$\n\n");
  printf("<!--SUMMARY_END--></FONT></B>\n");
}

void Bp3Result(const string result)
{
  printf("<B><FONT COLOR=\"#FF0000\"><!--SUMMARY_BEGIN-->\n\n");
  printf("$TEXT:Result: $$ A final result $$\n");
  printf("%s\n",result.c_str());
  printf("$$\n\n");
    printf("<!--SUMMARY_END--></FONT></B>\n\n");
}

// mathematical functions

double i0e(double x, int prec1, int prec2)
{
  /*
    Cephes Math Library Release 2.8:  June, 2000
    Copyright 1984, 1987, 2000 by Stephen L. Moshier
  */

  static double A[] =
    {
      -4.41534164647933937950E-18,
      3.33079451882223809783E-17,
      -2.43127984654795469359E-16,
      1.71539128555513303061E-15,
      -1.16853328779934516808E-14,
      7.67618549860493561688E-14,
      -4.85644678311192946090E-13,
      2.95505266312963983461E-12,
      -1.72682629144155570723E-11,
      9.67580903537323691224E-11,
      -5.18979560163526290666E-10,
      2.65982372468238665035E-9,
      -1.30002500998624804212E-8,
      6.04699502254191894932E-8,
      -2.67079385394061173391E-7,
      1.11738753912010371815E-6,
      -4.41673835845875056359E-6,
      1.64484480707288970893E-5,
      -5.75419501008210370398E-5,
      1.88502885095841655729E-4,
      -5.76375574538582365885E-4,
      1.63947561694133579842E-3,
      -4.32430999505057594430E-3,
      1.05464603945949983183E-2,
      -2.37374148058994688156E-2,
      4.93052842396707084878E-2,
      -9.49010970480476444210E-2,
      1.71620901522208775349E-1,
      -3.04682672343198398683E-1,
      6.76795274409476084995E-1
    };
  static double B[] =
    {
      -7.23318048787475395456E-18,
      -4.83050448594418207126E-18,
      4.46562142029675999901E-17,
      3.46122286769746109310E-17,
      -2.82762398051658348494E-16,
      -3.42548561967721913462E-16,
      1.77256013305652638360E-15,
      3.81168066935262242075E-15,
      -9.55484669882830764870E-15,
      -4.15056934728722208663E-14,
      1.54008621752140982691E-14,
      3.85277838274214270114E-13,
      7.18012445138366623367E-13,
      -1.79417853150680611778E-12,
      -1.32158118404477131188E-11,
      -3.14991652796324136454E-11,
      1.18891471078464383424E-11,
      4.94060238822496958910E-10,
      3.39623202570838634515E-9,
      2.26666899049817806459E-8,
      2.04891858946906374183E-7,
      2.89137052083475648297E-6,
      6.88975834691682398426E-5,
      3.36911647825569408990E-3,
      8.04490411014108831608E-1
    };

  double y;

  if (fabs(x) < DSMALL)
    return ONE;
  
  if( x < 0 )
    x = -x;
  if( x <= 8.0 )
  {
    y = (x/2.0) - 2.0;
    return( chbevl( y, A, 30, prec1 ) );
  }

  return(  chbevl( 32.0/x - 2.0, B, 25, prec2 ) / sqrt(x) );
}

double i1e( double x, int prec1, int prec2 )
{
  /*
    Cephes Math Library Release 2.8:  June, 2000 	
    Copyright 1985, 1987, 2000 by Stephen L. Moshier
  */

  static double A[] =
    {
      2.77791411276104639959E-18,
      -2.11142121435816608115E-17,
      1.55363195773620046921E-16,
      -1.10559694773538630805E-15,
      7.60068429473540693410E-15,
      -5.04218550472791168711E-14,
      3.22379336594557470981E-13,
      -1.98397439776494371520E-12,
      1.17361862988909016308E-11,
      -6.66348972350202774223E-11,
      3.62559028155211703701E-10,
      -1.88724975172282928790E-9,
      9.38153738649577178388E-9,
      -4.44505912879632808065E-8,
      2.00329475355213526229E-7,
      -8.56872026469545474066E-7,
      3.47025130813767847674E-6,
      -1.32731636560394358279E-5,
      4.78156510755005422638E-5,
      -1.61760815825896745588E-4,
      5.12285956168575772895E-4,
      -1.51357245063125314899E-3,
      4.15642294431288815669E-3,
      -1.05640848946261981558E-2,
      2.47264490306265168283E-2,
      -5.29459812080949914269E-2,
      1.02643658689847095384E-1,
      -1.76416518357834055153E-1,
      2.52587186443633654823E-1
    };

  static double B[] =
    {
      7.51729631084210481353E-18,
      4.41434832307170791151E-18,
      -4.65030536848935832153E-17,
      -3.20952592199342395980E-17,
      2.96262899764595013876E-16,
      3.30820231092092828324E-16,
      -1.88035477551078244854E-15,
      -3.81440307243700780478E-15,
      1.04202769841288027642E-14,
      4.27244001671195135429E-14,
      -2.10154184277266431302E-14,
      -4.08355111109219731823E-13,
      -7.19855177624590851209E-13,
      2.03562854414708950722E-12,
      1.41258074366137813316E-11,
      3.25260358301548823856E-11,
      -1.89749581235054123450E-11,
      -5.58974346219658380687E-10,
      -3.83538038596423702205E-9,
      -2.63146884688951950684E-8,
      -2.51223623787020892529E-7,
      -3.88256480887769039346E-6,
      -1.10588938762623716291E-4,
      -9.76109749136146840777E-3,
      7.78576235018280120474E-1
    };

  double y, z;

  z = fabs(x);
  if( z <= 8.0 )
  {
    y = (z/2.0) - 2.0;
    z = chbevl( y, A, 29, prec1 ) * z;
  }
  else
  {
    z = chbevl( 32.0/z - 2.0, B, 25, prec2 ) / sqrt(z);
  }
  if( x < 0.0 )
    z = -z;
  return( z );
}

double chbevl( double x, void* array, int n, int prec )
{
  /*
    Cephes Math Library Release 2.8:  June, 2000 	
    Copyright 1984, 1987, 2000 by Stephen L. Moshier
  */

  double b0, b1, b2, *p;
  int i;

  p = (double*) array + prec;
  b0 = *p++;
  b1 = 0.0;
  i = n - prec - 1;

  do
  {
    b2 = b1;
    b1 = b0;
    b0 = x * b1  -  b2  + *p++;
  }
  while( --i );

  return( 0.5*(b0-b2) );
}

double sim(const double x, int prec1, int prec2)
{
  return i1e(x,prec1,prec2)/i0e(x,prec1,prec2);
}

double simoverx(const double x, int prec1, int prec2)
{
  if (fabs(x) < DSMALL)
    return HALF;
  else
    return i1e(x,prec1,prec2)/(i0e(x,prec1,prec2)*x);
}
