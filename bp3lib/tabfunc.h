// tabulated sin,cos; I0,I1 on <0,200000> and exp on <-50,50>
// the (relative) precision in worst case is 10^-6 - don't use if highly precise results are requested
// approx. 2-15 x quicker than standard C implementations (depending on platform, compiler, optimization, function etc)
// Pavol Skubak, 2003-2007
// GNU GPL v2 licence

#ifndef TABFUNC_H
#define TABFUNC_H      1

#include <math.h>
#include "misc.h"
#include <fstream>
#include <iostream>
#include "tabfunctable.h"


using namespace std;


#if !defined PI_CONST
#define PI_CONST 1;
const double PIE = 3.14159265358979323846 ;
const double TWO_PI = 2*PIE ;
const double HALF_PI = 0.5*PIE ;
const double TWO_PI_INV = 1/TWO_PI ;
const double PI_INV = 1/PIE ;
#endif

template <typename realnum>
class TabFunc
{
  public:
    // ExpM(x) = exp(-x)
    realnum ExpM(realnum x)
    {
      if ( x<-50. || x>50. ) return exp(-x);
      return ExpM_nocheck(x);
    }

    // if you want to check whether the point is in <-50,50> region yourself or you now it must be in
    realnum ExpM_nocheck(realnum x)
    {
      poin = &t.Exparr[0] + (int) ( 250.5 + x*5. ); 
      x_diff = x - *poin;
      x_diff_sq = x_diff*x_diff;
      return ( *(poin+t.dimexp1)) * ( 1. - x_diff + 0.5*x_diff_sq -
          0.16666666666666*x_diff_sq*x_diff + 0.041666666666666*x_diff_sq*x_diff_sq ) ;
    }


    realnum Sin(realnum x)
    {
      ChargeTrig(x);
      return Sin_charged(x);
    }

    // use only if Cos(x) was called just before (saves time needed for table lookup)
    realnum Sin_charged(realnum x)
    {
      return ( *(poin+t.dimsincos1)) * ( 1. - 0.5*x_diff_sq + 0.041666666666666*x_diff_sq*x_diff_sq ) +
        ( *(poin+2*t.dimsincos1)) * ( x_diff - 0.16666666666666*x_diff_sq*x_diff ) ;
    }

    realnum Cos(realnum x)
    {
      ChargeTrig(x);
      return Cos_charged(x);
    }

    // use only if Sin(x) was called just before (saves time needed for table lookup)
    realnum Cos_charged(realnum x)
    {
      return (*(poin+2*t.dimsincos1)) * ( 1. - 0.5*x_diff_sq + 0.041666666666666*x_diff_sq*x_diff_sq )+
        (*(poin+t.dimsincos1)) * ( - x_diff + 0.16666666666666*x_diff_sq*x_diff ) ;
    }

    realnum I0e(realnum x)
    {
      if ( x < 0. ) return I0e_nocheck(-x);
      if ( x > 200000. ) return i0e(x,0,0);  
      return I0e_nocheck(x);
    }

    // if you want to check whether the point is in <0,50> region yourself or you now it must be in
    realnum I0e_nocheck(realnum x)
    {
      ChargeBess(x);
      return I0e_nocheck_charged(x);
    }

    // use only if I1e(x) was called just before (saves time needed for table lookup)
    realnum I0e_nocheck_charged(realnum x)
    {
      poin = &t.Bess0arr[i];
      return ( *(poin+t.dimbess1) + *(poin+2*t.dimbess1)*x_diff + *(poin+3*t.dimbess1)*x_diff_sq
          + *(poin+4*t.dimbess1)*x_diff_sq*x_diff + *(poin+5*t.dimbess1)*x_diff_sq*x_diff_sq );
    }

    realnum Sim(realnum x)
    {
      return I1e(x)/I0e(x);
    }

    realnum Simoverx(realnum x)
    {
      if (fabs(x) < 1.0e-7)
        return 0.5;
      else
        return I1e(x)/(x*I0e(x));
    }

    // use these two versions if you already computed I0e(x) and pass this value as the second argument
    realnum Sim_preI0(realnum x, realnum I0e_x_)
    {
      return I1e(x)/I0e_x_;
    }

    realnum Simoverx_preI0(realnum x,realnum I0e_x_)
    {
      if (fabs(x) < 1.0e-7)
        return 0.5;
      else
        return I1e(x)/(x*I0e_x_);
    }


    realnum I1e(realnum x)
    {
      if ( x<0. ) return -I1e_nocheck(-x);
      if ( x>200000. ) return i1e(x,0,0);
      return I1e_nocheck(x);
    }

    realnum I1e_nocheck(realnum x)
    {
      ChargeBess(x);
      return I1e_nocheck_charged(x);
    }

    // use only if I0e(x) was called just before (saves time needed for table lookup)
    realnum I1e_nocheck_charged(realnum x)
    {
      poin = &t.Bess1arr[i];
      return ( *(poin+t.dimbess1) + *(poin+2*t.dimbess1)*x_diff + *(poin+3*t.dimbess1)*x_diff_sq
          + *(poin+4*t.dimbess1)*x_diff_sq*x_diff + *(poin+5*t.dimbess1)*x_diff_sq*x_diff_sq );
    }

    realnum GetBess0coef(int i,int j)	{ return t.Bess0arr[i+j*t.dimbess1]; }
    realnum GetBess1coef(int i,int j)	{ return t.Bess1arr[i+j*t.dimbess1]; }
    realnum GetExpcoef(int i,int j)		{ return t.Exparr[i+j*t.dimexp1]; }
    realnum GetTrigcoef(int i,int j)	{ return t.SinCosarr[i+j*t.dimsincos1]; }
    unsigned GetBessdim(int i)			{ if (i==1) return t.dimbess1; else if (i==2) return t.dimbess2; else Error();   return 0; }
    unsigned GetExpdim(int i)			{ if (i==1) return t.dimexp1; else if (i==2) return t.dimexp2; else Error();   return 0; }
    unsigned GetTrigdim(int i)			{ if (i==1) return t.dimsincos1; else if (i==2) return t.dimsincos2; else Error();   return 0; }

  private:

    unsigned i;

    realnum x_diff, x_diff_sq;
    realnum *poin;
    static TabFuncTable<realnum> t;


    void ChargeBess(realnum x)
    {
      if ( x < 700. )
        if ( x < 30. )
          if ( x < 0. ) { x = 0.; i = 0; }
          else if ( x < 6. ) i = (int) ( .5 + x*5. );
          else i = (int) ( 30.5 + (x-6.)*1.25 );
        else
          if ( x < 140. ) i = (int) ( 60.5 + (x-30.)*0.25 );
          else i = (int) ( 88.5 + (x-140.)*0.05 );
      else
        if ( x < 30000. )
          if ( x < 4000. ) i = (int) ( 116.5 + (x-700.)*0.01 );
          else i = (int) ( 149.5 + (x-4000.)*0.002 );
        else
          if ( x < 86000. ) i = (int) ( 201.5 + (x-30000.)*0.001 );
          else i = (int) ( 257.5 + (x-86000.)*0.0005 );	  
      // using the fact that bessel0 and bessel1 are sampled in the same points here
      x_diff = x - t.Bess0arr[i];
      x_diff_sq = x_diff*x_diff;
    }

    void ChargeTrig(realnum x)
    {
      if ( x > TWO_PI )
        x -= TWO_PI * (realnum) ((int) (x*TWO_PI_INV));
      if ( x < 0. )
        x -= TWO_PI * (realnum) ((int) (x*TWO_PI_INV)-1);
      poin = &t.SinCosarr[0] + (int) ( 0.5 + x*4.7755491881566 );
      x_diff = x - *poin;
      x_diff_sq = x_diff*x_diff;
    }

    void Error() { cerr << "Wrong dimension number requested (only 1 or 2 allowed)" << endl; exit(1); }
};



#endif
