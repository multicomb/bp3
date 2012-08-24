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

#include <complex>
#include <algorithm>
#if defined (__sgi) && defined (_MIPS_ISA)
extern "C" {
#include <ctype.h>
}
#else
#include <cctype>
#include <cmath>
#endif
#include "gauss.h"
#include "likelihood.h"
#include "matrix.h"
#include "ccp4_fortran.h"
#include "csymlib.h"

#include "../mytimer.h"
extern TimerT Tinverse, TinverseGold;

extern "C" void FORTRAN_CALL ( DSYEVD, dsyevd,
    (char*, char*, int*, double*, int *,double*, double*, int*, int*, int*, int*),
    (char*, char*, int*, double*, int *,double*, double*, int*, int*, int*, int*),
    (char*, char*, int*, double*, int *,double*, double*, int*, int*, int*, int*));

extern "C" void
FORTRAN_CALL( ZHEEVD, zheevd,
    (char*, char*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, double*, int*, int*, int*, int*),
    (char*, char*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, double*, int*, int*, int*, int*),
    (char*, char*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, double*, int*, int*, int*, int*));

extern "C" 
{
    // LU decomoposition of a general matrix
    void dpotrf_(char* U, int *N, double* A, int* lda, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dpotri_(char* U, int *N, double* A, int* lda, int* INFO);

}
  Likelihood::Likelihood(Model &model, Crystal &crystal)
: mdl(model), xtal(crystal)
{
  // defaults
  likelihood = ZERO;
  verbose    = 0;
  warn       = false;
  hand       = 1;
  allin      = interpolate = stats = nooutput = only = false;
  protocol   = "PHASING";

  // default size for Gaussian integration parameters
  cweight.resize(0);
  pweight.resize(0);
  afweight.resize(0);
  sadweight.resize(0);
  maxoccshift     = TEN;
  maxnumbshift    = TEN;
  maxluzzshift    = TEN;
  maxbfacshift    = TEN;
  maxuanoshift    = 0.25;
  maxfpshift      = 0.5;
  maxfppshift     = 0.5;
  maxkscaleshift  = 0.25;
  maxbscaleshift  = 5.0;
  maxxyzshift     = maxxshift = maxyshift    = maxzshift = 0.05;   // ***NSP - should make it in terms of Orthogonal...
  rejectionprobability                       = 0.000000000000000000001;
  nosmall         = false;
  updatesigmah    = true;
  onlyanocalc     = false;
  fill            = true;
  beta            = 1.0;
  outputhcalc     = false;
  mtzfb           = "FB"; mtzpb              = "PHIB";
  mtzfom          = "FOM";
  mtzhla          = "HLA"; mtzhlb            = "HLB";
  mtzhlc          = "HLC"; mtzhld            = "HLD";
  mtzf            = "FPHASED"; mtzsigf       = "SIGFPHASED";
  mtzfa           = "FA";    mtzsigfa        = "SIGFA";
  mtzea           = "EA";    mtzsigea        = "SIGEA";
  mtzalpha        = "ALPHA";
  mtzfdiff        = "FDIFF"; mtzpdiff        = "PDIFF";
  mtzfcomb        = "FWT";   mtzpcomb        = "PHWT";
}

void Likelihood::checkocc()
{
  if ( ( (mode                   != "PHASE") || (target == "UNCO") )  )
  {
    bool setocc(false);
    // checks it occupancies are too high at user input
    for (unsigned i               = 0; i < mdl.atom.size(); i++)
      if ( (mdl.atom[i].Getocc()  > 0.5) && (xtal.sf.size() != 1) )
        setocc                    = true;

    if (verbose                  && setocc)
      Bp3Warning("Bp3likelihood::checkocc","Occupancies are high at user input, lowering them");

    if (setocc)
      for (unsigned i             = 0; i < mdl.atom.size(); i++)
        if (mdl.atom[i].Getocc()  > 0.5) 
          mdl.atom[i].Setocc(0.4);    
  }
}

void Likelihood::checkfrac(const vector<bool> &frac)
{
  // transform to fractional coordinates if orthogonal were given
  for (unsigned c                     = 0; c < xtal.cell.size(); c++)
    if (!frac[c])
      for (unsigned s                 = 0; s < mdl.site.size(); s++)
        // get crystal for the particular site
        for (unsigned a               = 0; a < mdl.atom.size(); a++)
          if ( (mdl.atom[a].nsite    == s) && (mdl.atom[a].crystal  == c) )
          {
            vector<double> f(3, ZERO);
            for (unsigned i           = 0; i < f.size(); i++)
            {
              for (unsigned j         = 0; j < f.size(); j++)
                f[i]                 += xtal.cell[c].or2frac[i][j]*mdl.site[s].x[j];
              mdl.site[s].Setx(i,f[i]);
            }
          }

  // make sure atoms are -1 < x < 1

  for (unsigned i                     = 0; i < mdl.site.size(); i++)
    for (unsigned x                   = 0; x < 3; x++)
      do
      {
        if (mdl.site[i].Getx(x)      >= 0.99)
          mdl.site[i].Setx(x, mdl.site[i].Getx(x) - ONE);
        else if (mdl.site[i].Getx(x) <= -ONE)
          mdl.site[i].Setx(x, mdl.site[i].Getx(x) + ONE);
      } while ( (mdl.site[i].Getx(x) >= 0.99) || (mdl.site[i].Getx(x) <= - ONE) );

}

void Likelihood::Setmaxshifts()
{
  // set the max shifts in fractional coordinates to depend on the maximum resolution
  // of the data.

  // determine maximum resolution
  double reso(100.0);
  for (unsigned d   = 0; d < xtal.sf.size(); d++)
    reso            = std::min((double)xtal.sf[d].hires,reso);

  for (unsigned i   = 0; i < 3; i++)
  {
    double shift(ZERO);
    for (unsigned j = 0; j < 3; j++)
      shift        += xtal.cell[0].or2frac[i][j]*reso;
    if (i          == 0)                                
      Setmaxx(shift);
    else if (i     == 1)
      Setmaxy(shift);
    else
      Setmaxz(shift);
  }

  if (verbose)
    printf("The maximum allowable shift in x, y and z is (resp.) %5.3f, %5.3f, %5.3f\n",maxxshift,maxyshift,maxzshift);
}


void Likelihood::checkatoms()
{

  // if occupancies are low, do not refine
  for (unsigned i               = 0; i < mdl.atom.size(); i++)
    if (mdl.atom[i].Getocc()    < 0.05)
    {
      mdl.site[mdl.atom[i].nsite].Setrefinex(0,false).Setrefinex(1,false).Setrefinex(2,false);
      mdl.atom[i].Setrefineocc(false).Setrefinebfac(false);
      printf("Not refining parameters for Atom %u\n",i+1);
    }
  printf("\n");
}

void Likelihood::resize(const bool outputmtz)
{
  // Allocate memory for variables.

  if (target            != "UNCO")
    for (unsigned d      = 0; d < xtal.sf.size(); d++)
      xtal.sf[d].refined = xtal.sf[d].refinead = false;

  if (target            != "SAD")
    for (unsigned d      = 0; d < xtal.sf.size(); d++)
      xtal.sf[d].refinesd  = false;

  if (target            != "MLHL")
    for (unsigned d      = 0; d < xtal.sf.size()  ; d++)
      xtal.sf[d].refinedmod = false;


  if (protocol          == "PHASECOMB")  
  {
    dLddmod.resize(xtal.sf.size());
    for (unsigned d      = 0; d < xtal.sf.size(); d++)
    {
      dLddmod[d].resize(xtal.sf[d].nbins);
      for (unsigned s    = 0; s < xtal.sf[d].nbins; s++)
        dLddmod[d][s].resize(xtal.sf[d].sfmodel.size(), ZERO);
    }
  }

  if ( (target          == "UNCO") && (dLddluz.size() != xtal.sf.size()) )
  {      
    // Non-Isomorphism parameter(s) 
    dLddluz.resize(xtal.sf.size());

    // Error in heavy atom model
    dLdeluz.resize(xtal.sf.size());

    // Error in anomalous difference
    dLdadluz.resize(xtal.sf.size());

    for (unsigned d      = 0; d < xtal.sf.size(); d++)
    {
      dLddluz[d].resize(xtal.sf[d].nbins, ZERO);
      dLdeluz[d].resize(xtal.sf[d].nbins, ZERO);
      dLdadluz[d].resize(xtal.sf[d].nbins, ZERO);
    }
  }

  if ( ( (target         == "MAD") || (target == "SAD") ) && (dLdsdluz.size() != xtal.sf.size()) )
  {
    // Luzzati Errors 
    //  dLddluz.resize(xtal.sf.size());
    dLdsdluz.resize(xtal.sf.size());

    // Matrices
    unsigned matrixsize(0);
    if (target           == "SAD")
      matrixsize          = 2;
    else
      for (unsigned d     = 0; d < xtal.sf.size(); d++)
        if (xtal.sf[d].anomalous)
          matrixsize     += 2;
        else
          matrixsize++;

    covmodel.resize(xtal.sf[0].nbins);
    covinvmodel.resize(xtal.sf[0].nbins);
    for (unsigned s       = 0; s < xtal.sf[0].nbins; s++)
    {
      covmodel[s].resize(matrixsize);
      covinvmodel[s].resize(matrixsize);
    }

    for (unsigned d       = 0; d < xtal.sf.size(); d++)
    {
      //      dLddluz[d].resize(xtal.sf[d].nbins, ZERO);
      dLdsdluz[d].resize(xtal.sf[d].nbins, ZERO);
    }
    unsigned row(2*matrixsize);
    recov.resize(row), recovinv.resize(row);
    redAdsd.resize(row), redAdsigh.resize(row);
    redAdsfpp.resize(row), redAddmod.resize(row);

    // inverse/eigenvalue filtering
    evalues.resize(row);

    lwork                 = 1 + 6*row + 2*row*row;
    liwork                = 3 + 5*row;
    lawork.resize(lwork);
    iwork.resize(liwork);    
  }

  // Pavol's functions

  if ( (target           == "MSRS") || (target == "PSAD") || (target == "PSD2") )
  {    
    dLdpluz.resize(xtal.pluz.size());

    for (unsigned p       = 0; p < xtal.pluz.size(); p++)
      xtal.refinep[p]     = (p < 2);

    if (!xtal.userpluz[1])
      xtal.refinep[1]     = true;

    if (protocol         == "PHASING")
    {
      if (refinexyz && refineocc && refinebfac)
        xtal.refinep[2]   = true;
      else
        xtal.refinep[2]   = false;    
    }
    else
      xtal.refinep[3]     = false;

    if (target           == "PSD2")
    {
      xtal.refinep[3]     = true;
      xtal.refinep[5]     = true;
    }


    unsigned ccd          = (target == "MSRS") ? 1 : 0;

    for (unsigned p       = 0; p < xtal.pluz.size(); p++)
    {
      if (!xtal.userpluz[p])
        xtal.pluz[p].resize(xtal.sf[ccd].nbins,xtal.defaultpluz[p]);

      dLdpluz[p].resize(xtal.sf[ccd].nbins,ZERO);
    }
  }

  // scale parameters
  dLdkscale.resize(xtal.sf.size(), ZERO);
  dLdbscale.resize(xtal.sf.size(), ZERO);

  // Real and imaginary F+ heavy atom structure factor (or just F if F+=F=)

  if (dLdAp.size()       != xtal.sf.size() )
  {
    dLdAp.resize(xtal.sf.size());
    dLdBp.resize(xtal.sf.size());

    for (unsigned d       = 0; d < xtal.sf.size(); d++)
    {
      dLdAp[d].resize(xtal.maxref, ZERO);
      dLdBp[d].resize(xtal.maxref, ZERO);
    }  

    // Real and imaginary F- heavy atom structure factor.
    // Memory is only allocated if anomalous scattering is present

    dLdAm.resize(xtal.sf.size());
    dLdBm.resize(xtal.sf.size());

    for (unsigned d       = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].anomalous)
      {	
        dLdAm[d].resize(xtal.maxref, ZERO);
        dLdBm[d].resize(xtal.maxref, ZERO);	
      }

    // Wilson's sigmaH parameter
    sumfpp.resize(xtal.sf.size());
    sumfpfpp.resize(xtal.sf.size());
    difffpp.resize(xtal.sf.size()), difffp.resize(xtal.sf.size());
    dLdsigmah.resize(xtal.sf.size());  
    dLdsumfpp.resize(xtal.sf.size());
    dLdsumfpfpp.resize(xtal.sf.size());

    for (unsigned d       = 0; d < xtal.sf.size(); d++)
    {
      sumfpp[d].resize(xtal.sf.size());
      sumfpfpp[d].resize(xtal.sf.size());
      difffpp[d].resize(xtal.sf.size());
      difffp[d].resize(xtal.sf.size());

      dLdsigmah[d].resize(xtal.sf[d].nbins, ZERO);
      dLdsumfpp[d].resize(xtal.sf[d].nbins, ZERO);
      dLdsumfpfpp[d].resize(xtal.sf[d].nbins, ZERO);

      for (unsigned d1    = 0; d1 < xtal.sf.size(); d1++)
      {
        sumfpp[d][d1].resize(xtal.sf[d].nbins, ZERO);
        sumfpfpp[d][d1].resize(xtal.sf[d].nbins, ZERO);
        difffpp[d][d1].resize(xtal.sf[d].nbins, ZERO);
        difffp[d][d1].resize(xtal.sf[d].nbins, ZERO);
      }
    }

    // Atomic parameters
    dLdx.resize(mdl.site.size());

    for (unsigned s       = 0; s < mdl.site.size(); s++)
      dLdx[s].resize(3, ZERO);  

    dLdocc.resize(mdl.atom.size(), ZERO);
    dLdnumb.resize(mdl.atom.size(), ZERO);

    dLdbiso.resize(mdl.atom.size(), ZERO);
    dLduaniso.resize(mdl.atom.size());

    for (unsigned a       = 0; a < mdl.atom.size(); a++)
      dLduaniso[a].resize(6, ZERO);

    dLdfp.resize(mdl.form.size());
    dLdfpp.resize(mdl.form.size());

    for (unsigned f       = 0; f < mdl.form.size(); f++)
    {
      dLdfp[f].resize(mdl.form[f].nwave, ZERO);
      dLdfpp[f].resize(mdl.form[f].nwave, ZERO);
    }

    // variables needed in structure factor calculation
    hrot.resize(xtal.sg.NSYMP);
    htran.resize(xtal.sg.NSYMP);
    for (unsigned s       = 0; s < hrot.size(); s++)
      hrot[s].resize(3);
  }

  // fom
  afom.resize(xtal.sf[0].nbins);
  cfom.resize(xtal.sf[0].nbins);

  anshl.resize(xtal.sf.size());
  cnshl.resize(xtal.sf.size());
  for (unsigned d        = 0; d < xtal.sf.size(); d++)
  {
    anshl[d].resize(xtal.sf[d].nbins, 0);
    cnshl[d].resize(xtal.sf[d].nbins, 0);
  }

  if (mode             != "EVALUES")
  {

    if ( (protocol       != "PHASECOMB") )
      if (outputmtz)
        xtal.Setselected("PHASE");
      else
        xtal.Setselected();

    // ***NSP
    if ( (protocol       == "PHASING") ||
        (protocol       == "PHASECOMB") && (target != "MLHL") )
      if (!updatesigmah)
        directsfcalc(false,!outputmtz);
      else
        directsfcalc(false,true);
  }

  xtal.Setnormalization().SetLuzzati();

  xtal.print();
}

void Likelihood::Settarget(const string targin)
{
  if ( (targin == "UNCO")    ||
      (targin == "SAD")     ||
      (targin == "MAD")     ||
      (targin == "MSRS")    ||
      (targin == "MLHL")    ||
      (targin == "PSD2")    ||
      (targin == "PSAD")      )   
    target      = targin;
  else
    Bp3Error("Crystal::Settarget", "Unknown target function");
}

double Likelihood::fp(const unsigned a, const unsigned w, const unsigned r)
{
  // give fp for a given atom, reflection and wavelength

  unsigned c(mdl.atom[a].crystal), f(mdl.atom[a].nform);
  double stolsq(xtal.cell[c].stolsq[r]);

  double normalscat((double)(mdl.form[f].a[0]*
        tab.ExpM(mdl.form[f].b[0]*stolsq) +
        mdl.form[f].a[1]*
        tab.ExpM(mdl.form[f].b[1]*stolsq) +
        mdl.form[f].a[2]*
        tab.ExpM(mdl.form[f].b[2]*stolsq) +
        mdl.form[f].a[3]*
        tab.ExpM(mdl.form[f].b[3]*stolsq) +
        mdl.form[f].c));

  return (normalscat + mdl.form[f].fp[w]);
}

template <class T> void copyzero(vector <T> &in, 
    vector <T> &out)
{
  out.resize(in.size(),0);
}

template <class T> void copyzero(vector< vector <T> > &in, 
    vector< vector <T> > &out)
{
  out.resize(in.size());
  for (unsigned int i = 0; i < in.size(); i++)
    copyzero(in[i],out[i]);
}

template <class T> void copyzero(vector<vector<vector <T> > > &in,
    vector<vector<vector <T> > > &out)
{
  out.resize(in.size());
  for (unsigned int i=0; i < in.size(); i++)
    copyzero(in[i],out[i]);
}

template <class T> void vectoradd(vector <T> &a,
    vector <T> &b)
{
  for (unsigned int i=0; i < a.size(); i++)
    a[i] += b[i];
}
template <class T> void vectoradd(vector<vector <T> > &a,
    vector<vector <T> > &b)
{
  for (unsigned int i=0; i < a.size(); i++)
    vectoradd(a[i],b[i]);
}

template <class T> void vectoradd(vector<vector<vector <T> > > &a,
    vector<vector<vector <T> > > &b)
{
  for (unsigned int i=0; i < a.size(); i++)
    vectoradd(a[i],b[i]);
}

template <class T> void vectorprint(const char *s, vector<T> &a)
{
  ofstream f;
  f.open(s,ios::app);
  f << s << endl;
  for (unsigned int i=0; i<a.size(); i++)
    f << a[i] << endl;
  f.close();
}

template <class T> void vectorprint(const char *s, vector<vector<T> > &a)
{
  for (unsigned int i=0; i<a.size(); i++)
    vectorprint(s,a[i]);
}

template <class T> void vectorprint(const char *s, vector<vector<vector<T> > > &a)
{
  for (unsigned int i=0; i<a.size(); i++)
    vectorprint(s,a[i]);
}

void Likelihood::directsfcalc(const bool deriv, const bool calcsigmah)
{
  // Calculates structure factors and first derivatives of a TARGET function
  // wrt parameters via direct summation.

  // The subroutine needs the derivative of the TARGET function
  // wrt the calculated structure amplitude and pcalc and updated
  // structure factor amplitude and pcalcs to calculate the required
  // derivatives.

  // The atomic coordinates are assumed to be in Fractional coordinates.
  // only sum over primitive symmetry operations
  double nmult((double)(xtal.sg.NSYM/xtal.sg.NSYMP));

  if (deriv)
  {
    // initialize first order partial derivatives
    for (unsigned a                        = 0; a < mdl.atom.size(); a++)
    {
      dLdocc[a]                            = ZERO;
      dLdbiso[a]                           = ZERO;

      for (unsigned b                      = 0; b < 6; b++)
        dLduaniso[a][b]                    = ZERO;
    }

    for (unsigned f                        = 0; f < mdl.form.size(); f++)
      for (unsigned w                      = 0; w < mdl.form[f].nwave; w++)
      {
        dLdfp[f][w]                        = ZERO;
        dLdfpp[f][w]                       = ZERO;
      }

    for (unsigned s                        = 0; s < mdl.site.size(); s++)
      for (unsigned k                      = 0; k < 3; k++)
        dLdx[s][k]                         = ZERO;
  }
  else if (calcsigmah)
    for (unsigned d                        = 0; d < xtal.sf.size(); d++)
      for (unsigned s                      = 0; s < xtal.sf[d].nbins; s++)
      {
        xtal.sf[d].sigmah[s]               = ZERO;

        for (unsigned d2                 = 0; d2 < sumfpp.size(); d2++)
        {
          sumfpp[d][d2][s]               = ZERO;
          difffpp[d][d2][s]              = ZERO;
          difffp[d][d2][s]               = ZERO;
        }
      }


#pragma omp parallel default(none) \
  shared(nmult)
  {
    vector<vector<work> >   omp_hrot(hrot); // need a local copy for omp
    // work defined in misc
    vector<work>            omp_htran(htran);
    vector<double>          omp_dLdocc;
    vector<double>          omp_dLdbiso;
    vector<vector<double> > omp_dLdx;
    vector<vector<double> > omp_dLduaniso(dLduaniso);

    copyzero(dLdocc,        omp_dLdocc);
    copyzero(dLdbiso,       omp_dLdbiso);
    copyzero(dLdx,          omp_dLdx);
    copyzero(dLduaniso,     omp_dLduaniso);

    //xtal.sf[d].sigmah[sh]     += scat*scat*btimeso*btimeso;
    //vector<Sfdata> omp_sf;

    vector<vector<double> > omp_sf_sigmah;

    vector< vector< vector<double> > > omp_sumfpp;
    vector< vector< vector<double> > > omp_difffpp;
    vector< vector< vector<double> > > omp_difffp;

    if (calcsigmah)
    {

      omp_sf_sigmah.resize(xtal.sf.size());
      for (unsigned int i=0; i < xtal.sf.size(); i++)
        omp_sf_sigmah[i].resize(xtal.sf[i].nbins,ZERO);

      copyzero(sumfpp,omp_sumfpp);
      copyzero(difffpp,omp_difffpp);
      copyzero(difffp,omp_difffp);
#if 0
      omp_sumfpp.resize(sumfpp.size());
      for (unsigned i=0; i < omp_sumfpp.size(); i++)
      {
        omp_sumfpp[i].resize(sumfpp[i].size());
        for (unsigned int j=0; j<omp_sumfpp[i].size(); j++)
          omp_sumfpp[i][j].resize(sumfpp[i][j].size(),ZERO);
      }

      omp_difffpp.resize(difffpp.size());
      for (unsigned i=0; i < omp_difffpp.size(); i++)
      {
        omp_difffpp[i].resize(difffpp[i].size());
        for (unsigned int j=0; j<omp_difffpp[i].size(); j++)
          omp_difffpp[i][j].resize(difffpp[i][j].size(),ZERO);
      }

      omp_difffp.resize(difffp.size());
      for (unsigned i=0; i < omp_difffp.size(); i++)
      {
        omp_difffp[i].resize(difffp[i].size());
        for (unsigned int j=0; j<omp_difffp[i].size(); j++)
          omp_difffp[i][j].resize(difffp[i][j].size(),ZERO);
      }
#endif

    }

    TabFunc<double> tab;

    vector<double> acalc(xtal.sf.size());
    vector<double> bcalc(xtal.sf.size());
    vector<double> anoacalc(xtal.sf.size());
    vector<double> anobcalc(xtal.sf.size());

#pragma omp for schedule(dynamic) nowait
    for (int r                          = 0; r < (int)xtal.maxselref; r++)
    {
      if (xtal.selected[r])
      {

        for (unsigned d                      = 0; d < xtal.sf.size(); d++)
        {
          acalc[d]                           = ZERO;
          bcalc[d]                           = ZERO;
          anoacalc[d]                        = ZERO;
          anobcalc[d]                        = ZERO;
        }

        double h(xtal.miller[r][0]), k(xtal.miller[r][1]), l(xtal.miller[r][2]);

        for (unsigned sym                    = 0; sym < xtal.sg.NSYMP; sym++)
        {
          omp_hrot[sym][0]                    = ((double)(xtal.sg.symrot[sym][0][0])*h +
              (double)(xtal.sg.symrot[sym][0][1])*k +
              (double)(xtal.sg.symrot[sym][0][2])*l);

          omp_hrot[sym][1]                    = ((double)(xtal.sg.symrot[sym][1][0])*h +
              (double)(xtal.sg.symrot[sym][1][1])*k +
              (double)(xtal.sg.symrot[sym][1][2])*l);

          omp_hrot[sym][2]                    = ((double)(xtal.sg.symrot[sym][2][0])*h +
              (double)(xtal.sg.symrot[sym][2][1])*k +
              (double)(xtal.sg.symrot[sym][2][2])*l);

          omp_htran[sym]                      = (xtal.sg.symtran[sym][0]*h +
              xtal.sg.symtran[sym][1]*k +
              xtal.sg.symtran[sym][2]*l);
        }

        for (unsigned a                      = 0; a < mdl.atom.size(); a++)
        {
          double costotal(ZERO), sintotal(ZERO);

          double dcosdx[3]                   = {ZERO, ZERO, ZERO};
          double dsindx[3]                   = {ZERO, ZERO, ZERO};

          // Atom's site number
          unsigned s(mdl.atom[a].nsite);

          for (unsigned sym                  = 0; sym < xtal.sg.NSYMP; sym++)
          {
            // Generate all symmetry elements
            double arg((omp_hrot[sym][0]*mdl.site[s].x[0] +
                  omp_hrot[sym][1]*mdl.site[s].x[1] +
                  omp_hrot[sym][2]*mdl.site[s].x[2] +
                  omp_htran[sym])*TWOPI);

            double cosarg(tab.Cos(arg));
            double sinarg(tab.Sin_charged(arg));

            // wwvv checked
            costotal                        += cosarg;
            sintotal                        += sinarg;

            if (deriv)
            {
              // wwvv checked
              dcosdx[0]                     -= sinarg*omp_hrot[sym][0];
              dsindx[0]                     += cosarg*omp_hrot[sym][0];
              dcosdx[1]                     -= sinarg*omp_hrot[sym][1];
              dsindx[1]                     += cosarg*omp_hrot[sym][1];
              dcosdx[2]                     -= sinarg*omp_hrot[sym][2];
              dsindx[2]                     += cosarg*omp_hrot[sym][2];
            }
          }

          // the crystal the atom belongs to
          unsigned c(mdl.atom[a].crystal);

          double stolsq(xtal.cell[c].stolsq[r]);

          double bfactor;

          double coefbano[6]                 = {ZERO, ZERO, ZERO, ZERO, ZERO, ZERO};

          if (mdl.atom[a].isotropic)
            bfactor                          = exp(-mdl.atom[a].biso*stolsq);
          else
          {
            coefbano[0]                      = h*h*xtal.cell[c].astar*xtal.cell[c].astar;
            coefbano[1]                      = TWO*h*k*xtal.cell[c].astar*xtal.cell[c].bstar;
            coefbano[2]                      = TWO*h*l*xtal.cell[c].astar*xtal.cell[c].cstar;
            coefbano[3]                      = k*k*xtal.cell[c].bstar*xtal.cell[c].bstar;
            coefbano[4]                      = TWO*k*l*xtal.cell[c].bstar*xtal.cell[c].cstar;
            coefbano[5]                      = l*l*xtal.cell[c].cstar*xtal.cell[c].cstar;
            bfactor                          = tab.ExpM(TWOPI2*std::max(coefbano[0]*
                  mdl.atom[a].uaniso[0] +
                  coefbano[1]*
                  mdl.atom[a].uaniso[1]+
                  coefbano[2]*
                  mdl.atom[a].uaniso[2]+
                  coefbano[3]*
                  mdl.atom[a].uaniso[3]+
                  coefbano[4]*
                  mdl.atom[a].uaniso[4]+
                  coefbano[5]*
                  mdl.atom[a].uaniso[5],
                  ZERO));
          }

          double btimeso(bfactor*mdl.atom[a].occ);

          //wwvv checked
          costotal                          *= bfactor*nmult;
          sintotal                          *= bfactor*nmult;

          if (deriv)
          {
            double temp(TWOPI*btimeso*nmult);
            //wwvv checked
            dcosdx[0]                       *= temp;
            dsindx[0]                       *= temp;
            dcosdx[1]                       *= temp;
            dsindx[1]                       *= temp;
            dcosdx[2]                       *= temp;
            dsindx[2]                       *= temp;
          }

          // Calculate the normal scattering factor for atom and resolution

          unsigned f(mdl.atom[a].nform);

          double normalscat((double)(mdl.form[f].a[0]*
                tab.ExpM(mdl.form[f].b[0]*stolsq) +
                mdl.form[f].a[1]*
                tab.ExpM(mdl.form[f].b[1]*stolsq) +
                mdl.form[f].a[2]*
                tab.ExpM(mdl.form[f].b[2]*stolsq) +
                mdl.form[f].a[3]*
                tab.ExpM(mdl.form[f].b[3]*stolsq) +
                mdl.form[f].c));

          // calculate anomalous scattering for each wavelength
          for (unsigned w                    = 0; w < mdl.form[f].fp.size(); w++)
          {
            unsigned d(xtal.dataset[c][w]);

            if (xtal.sf[d].use(r))
            {
              unsigned sh(xtal.bin(d,r));
              double scat(normalscat + mdl.form[f].fp[w]);

              if (onlyanocalc)
                scat                         = ZERO;

              double dReFpdo(scat*costotal), dImFpdo(scat*sintotal);

              if (!deriv)
              {
                //wwvv checked
                acalc[d]                    += dReFpdo*mdl.atom[a].occ;
                bcalc[d]                    += dImFpdo*mdl.atom[a].occ;

                if (calcsigmah)
                  // wwvv problem  solved
                  //xtal.sf[d].sigmah[sh]     += scat*scat*btimeso*btimeso;
                  omp_sf_sigmah[d][sh]        += scat*scat*btimeso*btimeso;

                if (xtal.sf[d].anomalous)
                {
                  double fpp(mdl.form[f].fpp[w]);
                  //wwvv checked
                  anoacalc[d]               -= fpp*sintotal*mdl.atom[a].occ;
                  anobcalc[d]               += fpp*costotal*mdl.atom[a].occ;

                  if (calcsigmah)
                  {
                    double fp(mdl.form[f].fp[w]*btimeso);
                    //wwvv checked
                    fpp                     *= btimeso;
                    // wwvv problem solved
                    //xtal.sf[d].sigmah[sh]   += fpp*fpp;
                    omp_sf_sigmah[d][sh]        += fpp*fpp;

                    // wwvv problem solved
                    if (sumfpp.size())
                    {
                      // wwvv problem solved
                      omp_sumfpp[d][d][sh]     += TWO*fpp*fpp;

                      for (unsigned w2      = w+1; w2 < mdl.form[f].fp.size(); w2++)
                      {
                        unsigned d2(xtal.dataset[c][w2]);
                        double fp2(mdl.form[f].fp[w2]*btimeso);
                        double fpp2(mdl.form[f].fpp[w2]*btimeso);
                        // wwvv problem solved  3*
                        omp_difffp[d][d2][sh]  += scat*(fp2 - fp)*btimeso;
                        omp_difffpp[d][d2][sh] += fpp*(fpp2 - fpp);
                        omp_sumfpp[d][d2][sh]  += fpp*(fpp2 + fpp);
                      }
                    }
                  }
                }
              }
              else
              {
                double dReFdx[3]             = {scat*dcosdx[0], scat*dcosdx[1], scat*dcosdx[2]};
                double dImFdx[3]             = {scat*dsindx[0], scat*dsindx[1], scat*dsindx[2]};

                if (xtal.sf[d].anomalous)
                {
                  //wwvv checked
                  dReFpdo                   -= mdl.form[f].fpp[w]*sintotal;
                  dImFpdo                   += mdl.form[f].fpp[w]*costotal;

                  for (unsigned i           = 0; i < 3; i++)
                  {
                    //wwvv checked
                    dReFdx[i]              -= mdl.form[f].fpp[w]*dsindx[i];
                    dImFdx[i]              += mdl.form[f].fpp[w]*dcosdx[i];
                  }
                }

                double mult(TWO*bfactor*btimeso*xtal.sg.NSYMP/

                    ((double)(xtal.sf[d].nshl[sh])));

                double temp(scat*scat*mult);

                // wwvv  problem solved: summation in dLdocc[a]
                omp_dLdocc[a]                += (dLdAp[d][r]*dReFpdo +
                    dLdBp[d][r]*dImFpdo +
                    dLdsigmah[d][sh]*temp);

                // wwvv problem solved with dLdx: 
                omp_dLdx[s][0]               += (dLdAp[d][r]*dReFdx[0] +
                    dLdBp[d][r]*dImFdx[0]);

                omp_dLdx[s][1]               += (dLdAp[d][r]*dReFdx[1] +
                    dLdBp[d][r]*dImFdx[1]);

                omp_dLdx[s][2]               += (dLdAp[d][r]*dReFdx[2] +
                    dLdBp[d][r]*dImFdx[2]);

                if (mdl.atom[a].isotropic)
                  // wwvv problem solved with dLdbiso
                  omp_dLdbiso[a]            -= mdl.atom[a].occ*stolsq*(dLdAp[d][r]*dReFpdo+
                      dLdBp[d][r]*dImFpdo+
                      dLdsigmah[d][sh]*
                      temp);
                else
                  for (unsigned u            = 0; u < 6; u++)
                    // wwvv problem solved with dLduaniso
                    omp_dLduaniso[a][u]      -= (TWOPI2*mdl.atom[a].occ*coefbano[u]*
                        (dLdAp[d][r]*dReFpdo +
                         dLdBp[d][r]*dImFpdo +
                         dLdsigmah[d][sh]*temp));

                if (xtal.sf[d].anomalous)
                {
                  double amult(TWO*bfactor*btimeso*
                      ((double)xtal.sg.NSYMP)/
                      ((double)(xtal.sf[d].anonshl[sh])));

                  double fpp(mdl.form[f].fpp[w]);
                  double dReFmdo(scat*costotal + fpp*sintotal);
                  double dImFmdo(scat*sintotal - fpp*costotal);

                  double atemp(fpp*fpp*mult);

                  // wwvv  problem solved: summation in dLdocc[a]
                  omp_dLdocc[a]              += (dLdAm[d][r]*dReFmdo +
                      dLdBm[d][r]*dImFmdo +
                      fpp*fpp*dLdsumfpp[d][sh]*amult +
                      dLdsigmah[d][sh]*atemp);

                  // wwvv problem solved with dLdx: 
                  omp_dLdx[s][0]             += (dLdAm[d][r]*(scat*dcosdx[0]+
                        fpp*dsindx[0])+
                      dLdBm[d][r]*(scat*dsindx[0]-
                        fpp*dcosdx[0]));

                  omp_dLdx[s][1]             += (dLdAm[d][r]*(scat*dcosdx[1]+
                        fpp*dsindx[1])+
                      dLdBm[d][r]*(scat*dsindx[1]-
                        fpp*dcosdx[1]));

                  omp_dLdx[s][2]             += (dLdAm[d][r]*(scat*dcosdx[2]+
                        fpp*dsindx[2])+
                      dLdBm[d][r]*(scat*dsindx[2]-
                        fpp*dcosdx[2]));

                  if (mdl.atom[a].isotropic)
                    // wwvv problem solved with dLdbiso
                    omp_dLdbiso[a]           -= (mdl.atom[a].occ*stolsq*
                        (dLdAm[d][r]*dReFmdo +
                         dLdBm[d][r]*dImFmdo +
                         fpp*fpp*dLdsumfpp[d][sh]*amult +
                         dLdsigmah[d][sh]*atemp));
                  else
                    for (unsigned u          = 0; u < 6; u++)
                      // wwvv problem solved with dLduaniso
                      omp_dLduaniso[a][u]    -= (TWOPI2*mdl.atom[a].occ*coefbano[u]*
                          (dLdAm[d][r]*dReFmdo +
                           dLdBm[d][r]*dImFmdo +
                           fpp*fpp*dLdsumfpp[d][sh]*amult +
                           dLdsigmah[d][sh]*atemp));
                }
              }
            }
          }
        } // end for(unsigned a=0; ...)

        if (!deriv) // store amplitude and phase
          for (unsigned d                    = 0; d < xtal.sf.size(); d++)
          {
            double atemp(acalc[d] + anoacalc[d]), btemp(bcalc[d] + anobcalc[d]);

            xtal.sf[d].fcalcp[r]             = sqrt(atemp*atemp + btemp*btemp);
            xtal.sf[d].pcalcp[r]             = atan2(btemp, atemp);

            if (xtal.sf[d].anomalous)
            {
              atemp                          = acalc[d] - anoacalc[d];
              btemp                          = bcalc[d] - anobcalc[d];
              xtal.sf[d].fcalcm[r]           = sqrt(atemp*atemp + btemp*btemp);
              xtal.sf[d].pcalcm[r]           = atan2(btemp, atemp);
            }
          }
      } //if xtal.selected...
    }  // end omp for(r=0;...)

#pragma omp critical
    {

      vectoradd(dLdocc,    omp_dLdocc);
      vectoradd(dLdbiso,   omp_dLdbiso);
      vectoradd(dLdx,      omp_dLdx);
      vectoradd(dLduaniso, omp_dLduaniso);

#if 0
      for (unsigned int i=0; i < dLdocc.size(); i++)
        dLdocc[i] += omp_dLdocc[i];

      for(unsigned int i=0; i < dLdbiso.size(); i++)
        dLdbiso[i] += omp_dLdbiso[i];

      for (unsigned int i=0; i < dLdx.size(); i++)
        for (unsigned int j=0; j < dLdx[i].size(); j++)
          dLdx[i][j] += omp_dLdx[i][j];

      for (unsigned int i=0; i < dLduaniso.size(); i++)
        for (unsigned int j=0; j < dLduaniso[i].size(); j++)
          dLduaniso[i][j] += omp_dLduaniso[i][j];
#endif

      if (calcsigmah)
      {
        for (unsigned int i=0; i<omp_sf_sigmah.size(); i++)
          for (unsigned int j=0; j<omp_sf_sigmah[i].size(); j++)
            xtal.sf[i].sigmah[j] += omp_sf_sigmah[i][j];

        vectoradd(difffp, omp_difffp);
        vectoradd(difffpp,omp_difffpp);
        vectoradd(sumfpp, omp_sumfpp);
#if 0
        for (unsigned i=0; i < difffp.size(); i++)
          for (unsigned j=0; j < difffp[i].size(); j++)
            for (unsigned k=0; k < difffp[i][j].size(); k++)
              difffp[i][j][k]  += omp_difffp[i][j][k];

        for (unsigned i=0; i < difffpp.size(); i++)
          for (unsigned j=0; j < difffpp[i].size(); j++)
            for (unsigned k=0; k < difffpp[i][j].size(); k++)
              difffpp[i][j][k] += omp_difffpp[i][j][k];

        for (unsigned i=0; i < omp_sumfpp.size(); i++)
          for (unsigned j=0; j < omp_sumfpp[i].size(); j++)
            for (unsigned k=0; k < omp_sumfpp[i][j].size(); k++)
              sumfpp[i][j][k]  += omp_sumfpp[i][j][k];
#endif
      }
    }
  }    // end omp parallel

  /*
     vectorprint("dLdocc.txt",dLdocc);
     vectorprint("dLdbiso.txt",dLdbiso);
     vectorprint("dLdx.txt",dLdx);
     vectorprint("dLduaniso.txt",dLduaniso);
     vectorprint("sumfpp.txt",sumfpp);
     vectorprint("difffpp.txt",difffpp);
     vectorprint("difffp.txt",difffp);
#if 0
{
fstream f;
f.open("xtalsfsigmah.txt",ios::app);
for (unsigned int i=0; i<xtal.sf.size(); i++)
for (unsigned int j=0; j<xtal.sf[i].sigmah.size(); j++)
f << xtal.sf[i].sigmah[j] << endl;
f.close();
}
#endif
*/

if (!deriv)
  for (unsigned d                        = 0; d < xtal.sf.size()   ; d++)
  for (unsigned s                      = 0; s < xtal.sf[d].nbins ; s++)
{
  if (xtal.sf[d].nshl[s]            && calcsigmah)
    xtal.sf[d].sigmah[s]            *= ((double)(xtal.sg.NSYMP)/
        ((double)(xtal.sf[d].nshl[s])));

  if (xtal.sf[d].anomalous          && calcsigmah &&
      ( (target                     == "SAD")  || (target == "MAD") ||
        (target                     == "PSAD") || (target == "MSRS") ||
        (target                     == "PSD2")) )
    if (xtal.sf[d].anonshl[s])
      for (unsigned d2               = 0; d2 < xtal.sf.size()  ; d2++)
      {
        sumfpp[d][d2][s]            *= ((double)(xtal.sg.NSYMP)/
            ((double)xtal.sf[d].anonshl[s]));
        difffpp[d][d2][s]           *= ((double)(xtal.sg.NSYMP)/
            ((double)xtal.sf[d].anonshl[s]));
        difffp[d][d2][s]            *= ((double)(xtal.sg.NSYMP)/
            ((double)xtal.sf[d].anonshl[s]));
      }
}
}

double Likelihood::calcanomcorrel() const
{
  const double FPPCUT(1.5);

  double resover                = 0.0;

  if (xtal.sf.size()            > 1)
    for (unsigned d1            = 0; d1 < xtal.sf.size(); d1++)
      for (unsigned d2          = d1+1; d2 < xtal.sf.size(); d2++)
        if (xtal.sf[d1].anomalous && xtal.sf[d2].anomalous &&
            (xtal.sf[d1].nxtal == xtal.sf[d2].nxtal) )
        {
          bool lowfpp(false);
          unsigned d            = ( (xtal.sf[d1].hires > xtal.sf[d2].hires) ?
              d1 : d2);

          vector<double> sumx(xtal.sf[d].nbins, ZERO);
          vector<double> sumx2(xtal.sf[d].nbins, ZERO);
          vector<double> sumy(xtal.sf[d].nbins, ZERO);
          vector<double> sumy2(xtal.sf[d].nbins, ZERO);
          vector<double> sumxy(xtal.sf[d].nbins, ZERO);
          vector<double> sum(xtal.sf[d].nbins, ZERO);

          for (unsigned f       = 0; f < mdl.form.size(); f++)
            for (unsigned w     = 0; w < mdl.form[f].nwave; w++)
            {
              unsigned c(mdl.form[f].Getcrystal());
              if ( (mdl.form[f].fpp[w] < FPPCUT) &&
                  ( (xtal.dataset[c][w] == d1) || (xtal.dataset[c][w] == d2) ) )
                lowfpp          = true;
            }


          for (unsigned r       = 0; r < xtal.maxref; r++)
            if (xtal.sf[d1].anouse(r) && xtal.sf[d2].anouse(r))
            {
              unsigned s(xtal.bin(d,r));

              sumx[s]          += xtal.sf[d1].dano(r);
              sumy[s]          += xtal.sf[d2].dano(r);
              sumxy[s]         += xtal.sf[d1].dano(r)*xtal.sf[d2].dano(r);
              sumx2[s]         += xtal.sf[d1].dano(r)*xtal.sf[d1].dano(r); 
              sumy2[s]         += xtal.sf[d2].dano(r)*xtal.sf[d2].dano(r);
              sum[s]           += ONE;
            }

          double totals(ZERO), totalcor(ZERO), correl(ZERO);
          printf("\n $TABLE: Signed anomalous difference correlation: %s and %s\n",xtal.sf[d1].name.c_str(),
              xtal.sf[d2].name.c_str());
          printf("$GRAPHS: Signed anomalous difference correlation vs. 1/Res^2 :N:4,6:\n");
          printf("$$\n Bin   HiRes  LoRes  1/Res^2   Refls  Correlation $$\n$$\n");
          double lastcorrel(ZERO), lastreso(ZERO);
          bool underthirty(false);
          for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
          {
            if (sum[s]          > ZERO)
              correl            = ((sum[s]*sumxy[s] - sumx[s]*sumy[s])/
                  (sqrt(sum[s]*sumx2[s] - sumx[s]*sumx[s])*
                   sqrt(sum[s]*sumy2[s] - sumy[s]*sumy[s])));
            else
              correl            = ZERO;

            double reso         = HALF/sqrt(xtal.sf[d].astolsq[s]);

            if ( (correl        < 0.300) && (lastcorrel > 0.300) && (!underthirty)
                && (!lowfpp) )
            {
              double tmp        = (lastreso-reso)/(lastcorrel-correl)*(0.3-correl) + reso;
              underthirty       = true;
              if ( (resover     > tmp) || (resover == ZERO) )
                resover         = tmp;
            }

            lastcorrel          = correl;
            lastreso            = reso;

            printf(" %2d  %6.2f  %6.2f  %7.5f  %5.0f  %7.5f\n",s+1,xtal.sf[d].Getlores(s),
                xtal.sf[d].Gethires(s), xtal.sf[d].astolsq[s],sum[s],correl);
            totalcor           += correl*sum[s];
            totals             += sum[s];
          }
          if (!underthirty     && (resover == ZERO) && !(lowfpp))
            if (correl         >= 0.3)
              resover           = xtal.sf[d].hires;
          printf("$$\nTOTAL\n");
          printf("     %6.2f  %6.2f           %5.0f  %7.5f\n\n",xtal.sf[d].Getlores(0),
              xtal.sf[d].Gethires(xtal.sf[d].nbins-1),totals,totalcor/totals);
        }
  return resover;
}

void Likelihood::shift(const double stepsize, const vector<double> &pars,
    const vector<double> &direction)
{
  int i(-1);

  // Apply direction and step size to the parameters

  for (unsigned d       = 0; d < xtal.sf.size()            ; d++)
  {
    if (xtal.sf[d].refinek && refinescale)
    {
      ++i; xtal.sf[d].Setkscale(pars[i] + xsshift(stepsize*direction[i], "KSCALE"));
    }
    if (xtal.sf[d].refineb && refinescale)
    {
      ++i; xtal.sf[d].Setbiso(pars[i] + xsshift(stepsize*direction[i], "BSCALE"));
    }

    if (xtal.sf[d].refined  && refineluzzati)
      for (unsigned s   = 0; s < xtal.sf[d].dluz.size()    ; s++)
      {
        ++i; xtal.sf[d].Setdluz(s,pars[i] + xsshift(stepsize*direction[i], "LUZZ"));
      }

    if (xtal.sf[d].refinee  && refineluzzati)
      for (unsigned s   = 0; s < xtal.sf[d].eluz.size()    ; s++)
      {
        ++i; xtal.sf[d].Seteluz(s,pars[i] + xsshift(stepsize*direction[i], "LUZZ"));
      }

    if (xtal.sf[d].refinead && refineluzzati)
      for (unsigned s   = 0; s < xtal.sf[d].adluz.size()   ; s++)
      {
        ++i; xtal.sf[d].Setadluz(s,pars[i] + xsshift(stepsize*direction[i], "LUZZ"));
      }

    if (xtal.sf[d].refinesd && refineluzzati)
      for (unsigned s   = 0; s < xtal.sf[d].sdluz.size()   ; s++)
      {
        ++i; xtal.sf[d].Setsdluz(s,pars[i] + xsshift(stepsize*direction[i], "LUZZ"));
      }

    if (xtal.sf[d].refinedmod && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].dmod.size()    ; s++)
        for (unsigned m   = 0; m < xtal.sf[d].dmod[s].size() ; m++)
        {
          ++i; xtal.sf[d].Setdmod(s, m, pars[i] + xsshift(stepsize*direction[i], "LUZZ"));
        }
  }

  for (unsigned p       = 0; p < xtal.pluz.size()   ; p++)
    if (xtal.refinep[p] && refineluzzati)
      for (unsigned s   = 0; s < xtal.pluz[p].size(); s++)
      {
        ++i; xtal.Setpluz(p,s,pars[i] + xsshift(stepsize*direction[i], "LUZZ"));
      }

  if (refinexyz)
    for (unsigned s     = 0; s < mdl.site.size()           ; s++)
      for (unsigned k   = 0; k < 3                         ; k++)
        if (mdl.site[s].Getrefinex(k))
        {
          ++i; mdl.site[s].Setx(k, pars[i] + xsshift(stepsize*direction[i], "XYZ"));
        }

  for (unsigned a       = 0; a < mdl.atom.size()           ; a++)
  {
    if (mdl.atom[a].Getrefinebfac() && refinebfac)
      if (mdl.atom[a].Getisotropic() )
      {
        ++i; mdl.atom[a].Setbiso(pars[i] + xsshift(stepsize*direction[i], "BFAC"));
      }
      else
        for (unsigned k = 0; k < 6                         ; k++)
        {
          ++i; mdl.atom[a].Setuaniso(pars[i] + xsshift(stepsize*direction[i], "UANO"), k);  
        }

    if (mdl.atom[a].Getrefineocc() && refineocc)
    {
      ++i; mdl.atom[a].Setocc(pars[i] + xsshift(stepsize*direction[i], "OCCU"));
    }

    if (mdl.atom[a].Getrefinenumb() && refinenumb)
    {
      ++i; mdl.atom[a].Setnumb(pars[i] + xsshift(stepsize*direction[i], "NUMB"));
    }
  }

  for (unsigned f       = 0; f < mdl.form.size()           ; f++)
    for (unsigned w     = 0; w < mdl.form[f].Getnwave()    ; w++)
    {
      if (mdl.form[f].Getrefinefp(w) && refinefp)
      {
        ++i; mdl.form[f].Setfp(w, pars[i] + xsshift(stepsize*direction[i], "FP"));
      }

      if (mdl.form[f].Getrefinefpp(w) && refinefpp)
      {
        ++i; mdl.form[f].Setfpp(w, pars[i] + xsshift(stepsize*direction[i], "FPP"));
      }
    }
}

double Likelihood::xsshift(const double shift, const string &type) const
{
  // dampen shifts - if required

  if (type     == "OCCU")
    return ( (fabs(shift) > maxoccshift)    ? sign(shift)*maxoccshift    : shift);
  else if (type == "NUMB")
    return ( (fabs(shift) > maxnumbshift)   ? sign(shift)*maxnumbshift   : shift);
  else if (type == "LUZZ")
    return ( (fabs(shift) > maxluzzshift)   ? sign(shift)*maxluzzshift   : shift);
  else if (type == "BFAC")
    return ( (fabs(shift) > maxbfacshift)   ? sign(shift)*maxbfacshift   : shift);
  else if (type == "XYZ")
    return ( (fabs(shift) > maxxyzshift)    ? sign(shift)*maxxyzshift    : shift);
  else if (type == "X")
    return ( (fabs(shift) > maxxshift)      ? sign(shift)*maxxshift      : shift);
  else if (type == "Y")
    return ( (fabs(shift) > maxyshift)      ? sign(shift)*maxyshift      : shift);
  else if (type == "Z")
    return ( (fabs(shift) > maxzshift)      ? sign(shift)*maxzshift      : shift);
  else if (type == "UANO")
    return ( (fabs(shift) > maxuanoshift)   ? sign(shift)*maxuanoshift   : shift);
  else if (type == "FP")
    return ( (fabs(shift) > maxfpshift)     ? sign(shift)*maxfpshift     : shift);
  else if (type == "FPP")
    return ( (fabs(shift) > maxfppshift)    ? sign(shift)*maxfppshift    : shift);
  else if (type == "KSCALE")
    return ( (fabs(shift) > maxkscaleshift) ? sign(shift)*maxkscaleshift : shift);
  else if (type == "BSCALE")
    return ( (fabs(shift) > maxbscaleshift) ? sign(shift)*maxbscaleshift : shift);
  else
    Bp3Error("Bp3likelihood::xsshift", "Type unknown");

  return ZERO;
}

vector<double> Likelihood::Getgradient() const
{
  // place gradients in a vector

  vector<double> grad;

  for (unsigned d           = 0; d < xtal.sf.size()          ; d++)
  {
    if (xtal.sf[d].refinek  && refinescale)
      grad.push_back(dLdkscale[d]);
    if (xtal.sf[d].refineb && refinescale)
      grad.push_back(dLdbscale[d]);

    if (xtal.sf[d].refined  && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].dluz.size()    ; s++)
        grad.push_back(dLddluz[d][s]);
    if (xtal.sf[d].refinee  && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].eluz.size()    ; s++)
        grad.push_back(dLdeluz[d][s]);
    if (xtal.sf[d].refinead && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].adluz.size()   ; s++)
        grad.push_back(dLdadluz[d][s]);
    if (xtal.sf[d].refinesd && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].sdluz.size()   ; s++)
        grad.push_back(dLdsdluz[d][s]);
    if (xtal.sf[d].refinedmod && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].dmod.size()    ; s++)
        for (unsigned m   = 0; m < xtal.sf[d].dmod[s].size() ; m++)
          grad.push_back(dLddmod[d][s][m]);
  }

  for (unsigned p         = 0; p < xtal.pluz.size()      ; p++)
    if (xtal.refinep[p] && refineluzzati)
      for (unsigned s     = 0; s < xtal.pluz[p].size()   ; s++)
        grad.push_back(dLdpluz[p][s]);

  if (refinexyz)
    for (unsigned s         = 0; s < mdl.site.size()       ; s++)
      for (unsigned k       = 0; k < 3                     ; k++)
        if (mdl.site[s].Getrefinex(k))
          grad.push_back(dLdx[s][k]);

  for (unsigned a         = 0; a < mdl.atom.size()       ; a++)
  {
    if (mdl.atom[a].Getrefinebfac() && refinebfac)
      if (mdl.atom[a].Getisotropic() )
        grad.push_back(dLdbiso[a]);
      else
      {
        // ***NSP - check to see if anisotropic matrix is positive definite 
        if (!mdl.atom[a].posdefu() )
          printf("Atom %u does not have a positive definite anisotropic matrix\n",a+1);

        for (int b        = 0; b < 6                     ; b++)
          grad.push_back(dLduaniso[a][b]);
      }

    if (mdl.atom[a].Getrefineocc() && refineocc)
      grad.push_back(dLdocc[a]);

    if (mdl.atom[a].Getrefinenumb() && refinenumb)
      grad.push_back(dLdnumb[a]);

  }

  for (unsigned f         = 0; f < mdl.form.size()       ; f++)
    for (unsigned w       = 0; w < mdl.form[f].Getnwave(); w++)
    {
      if (mdl.form[f].Getrefinefp(w)  && refinefp)
        grad.push_back(dLdfp[f][w]);
      if (mdl.form[f].Getrefinefpp(w) && refinefpp)
        grad.push_back(dLdfpp[f][w]);
    }

  return grad;
}

double Likelihood::Getstepmax(const vector<double> &direction) const
{
  double step(ZERO);
  for (unsigned i          = 0; i < direction.size(); i++)
    if (fabs(direction[i]) > EPSILON)
    {
      if (type[i]         == "LUZZ")
        step               = std::max(step, fabs(maxluzzshift/direction[i]));
      if (type[i]         == "XYZ")
        step               = std::max(step, fabs(maxxyzshift/direction[i]));
      /*
         if (type[i]         == "BFAC")
         step               = std::max(step, fabs(maxbfacshift/direction[i]));
         */
      if (type[i]         == "OCCU")
        step               = std::max(step, fabs(maxoccshift/direction[i]));
      if (type[i]         == "NUMB")
        step               = std::max(step, fabs(maxnumbshift/direction[i]));
      if (type[i]         == "KSCALE")
        step               = std::max(step, fabs(maxkscaleshift/direction[i]));
      if (type[i]         == "BSCALE")
        step               = std::max(step, fabs(maxbscaleshift/direction[i]));
    }

  step                     = between(ZERO,step,ONE);

  return step;
}

void Likelihood::Settype()
{
  type.resize(0);

  // place type of parameter in vector

  for (unsigned d           = 0; d < xtal.sf.size()          ; d++)
  {
    if (xtal.sf[d].refinek && refinescale)
      type.push_back("KSCALE");
    if (xtal.sf[d].refineb && refinescale)
      type.push_back("BSCALE");

    if (xtal.sf[d].refined && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].dluz.size()    ; s++)
        type.push_back("LUZZ");
    if (xtal.sf[d].refinee && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].eluz.size()    ; s++)
        type.push_back("LUZZ");
    if (xtal.sf[d].refinead && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].adluz.size()   ; s++)
        type.push_back("LUZZ"); 
    if (xtal.sf[d].refinesd && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].sdluz.size()   ; s++)
        type.push_back("LUZZ");

    if (xtal.sf[d].refinedmod && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].dmod.size()    ; s++)
        for (unsigned m   = 0; m < xtal.sf[d].dmod[s].size() ; m++)
          type.push_back("LUZZ");
  }

  for (unsigned p        = 0; p < xtal.pluz.size()           ; p++)
    if (xtal.refinep[p] && refineluzzati)
      for (unsigned s    = 0; s < xtal.pluz[p].size()        ; s++)
        type.push_back("LUZZ");

  if (refinexyz)
    for (unsigned s      = 0; s < mdl.site.size()            ; s++)
      for (unsigned k    = 0; k < 3                          ; k++)
        if (mdl.site[s].Getrefinex(k))
          type.push_back("XYZ");

  for (unsigned a      = 0; a < mdl.atom.size()            ; a++)
  {
    if (mdl.atom[a].Getrefinebfac() && refinebfac)
      if (mdl.atom[a].Getisotropic() )
        type.push_back("BFAC");
      else
      {
        // ***NSP - check to see if anisotropic matrix is positive definite 
        if (!mdl.atom[a].posdefu() )
        {
          printf("Atom %u", a+1);
          printf(" does not have a positive definite anisotropic matrix\n");
        }

        for (int b     = 0; b < 6                        ; b++)
          type.push_back("UANO");
      }

    if (mdl.atom[a].Getrefineocc() && refineocc)
      type.push_back("OCCU");

    if (mdl.atom[a].Getrefinenumb() && refinenumb)
      type.push_back("NUMB");
  }

  for (unsigned f      = 0; f < mdl.form.size()          ; f++)
    for (unsigned w    = 0; w < mdl.form[f].Getnwave()   ; w++)
    {
      if (mdl.form[f].Getrefinefp(w)  && refinefp)
        type.push_back("FP");
      if (mdl.form[f].Getrefinefpp(w) && refinefpp)
        type.push_back("FPP");
    }
}

vector<double> Likelihood::Setpars() const
{
  // Place original/refineable parameters in array pars

  vector<double> pars;

  for (unsigned d         = 0; d < xtal.sf.size()           ; d++)
  {
    if (xtal.sf[d].refinek && refinescale)
      pars.push_back(xtal.sf[d].Getkscale());

    if (xtal.sf[d].refineb && refinescale)
      pars.push_back(xtal.sf[d].Getbiso());

    if (xtal.sf[d].refined  && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].dluz.size()    ; s++)
        pars.push_back(xtal.sf[d].Getdluz(s));

    if (xtal.sf[d].refinee  && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].eluz.size()    ; s++)
        pars.push_back(xtal.sf[d].Geteluz(s));

    if (xtal.sf[d].refinead && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].adluz.size()   ; s++)
        pars.push_back(xtal.sf[d].Getadluz(s));

    if (xtal.sf[d].refinesd && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].sdluz.size()   ; s++)
        pars.push_back(xtal.sf[d].Getsdluz(s));

    if (xtal.sf[d].refinedmod && refineluzzati)
      for (unsigned s     = 0; s < xtal.sf[d].dmod.size()    ; s++)
        for (unsigned m   = 0; m < xtal.sf[d].dmod[s].size() ; m++)
          pars.push_back(xtal.sf[d].Getdmod(s,m));
  }

  for (unsigned p         = 0; p < xtal.pluz.size()          ; p++)
    if (xtal.refinep[p]  && refineluzzati)
      for (unsigned s     = 0; s < xtal.pluz[p].size()       ; s++)
        pars.push_back(xtal.Getpluz(p,s));

  if (refinexyz)
    for (unsigned s         = 0; s < mdl.site.size()         ; s++)
      for (unsigned k       = 0; k < 3                       ; k++)
        if (mdl.site[s].Getrefinex(k))
          pars.push_back(mdl.site[s].Getx(k));

  for (unsigned a         = 0; a < mdl.atom.size()         ; a++)
  {
    if (mdl.atom[a].Getrefinebfac() && refinebfac)
      if (mdl.atom[a].Getisotropic() )
        pars.push_back(mdl.atom[a].Getbiso());
      else
        for (unsigned b   = 0; b < 6                       ; b++)
          pars.push_back(mdl.atom[a].Getuaniso(b));

    if (mdl.atom[a].Getrefineocc() && refineocc)
      pars.push_back(mdl.atom[a].Getocc());      

    if (mdl.atom[a].Getrefinenumb() && refinenumb)
      pars.push_back(mdl.atom[a].Getnumb());      
  }

  for (unsigned f         = 0; f < mdl.form.size()         ; f++)
    for (unsigned w       = 0; w < mdl.form[f].Getnwave()  ; w++)
    {
      if (mdl.form[f].Getrefinefp(w)  && refinefp)
        pars.push_back(mdl.form[f].Getfp(w));
      if (mdl.form[f].Getrefinefpp(w) && refinefpp)
        pars.push_back(mdl.form[f].Getfpp(w));
    }

  return pars;
}

void Likelihood::writepdb() const
{
  for (unsigned c                      = 0; c < xtal.cell.size(); c++)
  {
    bool atomincrystal(false);
    for (unsigned a                    = 0; a < mdl.atom.size(); a++)
      if (c                           == mdl.atom[a].crystal)
        atomincrystal                  = true;

    if (atomincrystal)
    {
      char pdbnumber[3];

      if (c                            < 9)
        sprintf(pdbnumber,"%1u",c+1);
      else
        sprintf(pdbnumber,"%2u",c+1);

      string pdbname                   = output + "-" + pdbnumber + ".pdb";

      FILE *pdbfile                    = fopen(pdbname.c_str(),"w");

      fprintf(pdbfile,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %10s\n",
          xtal.cell[c].a, xtal.cell[c].b, xtal.cell[c].c,
          xtal.cell[c].alpha, xtal.cell[c].beta, xtal.cell[c].gamma,
          xtal.sg.Gethmname().c_str());
      fprintf(pdbfile,"SCALE1    %10.6f%10.6f%10.6f     %10.5f\n",
          xtal.cell[c].or2frac[0][0], xtal.cell[c].or2frac[0][1],
          xtal.cell[c].or2frac[0][2],0.0);
      fprintf(pdbfile,"SCALE2    %10.6f%10.6f%10.6f     %10.5f\n",
          xtal.cell[c].or2frac[1][0], xtal.cell[c].or2frac[1][1],
          xtal.cell[c].or2frac[1][2],0.0);
      fprintf(pdbfile,"SCALE3    %10.6f%10.6f%10.6f     %10.5f\n",
          xtal.cell[c].or2frac[2][0], xtal.cell[c].or2frac[2][1],
          xtal.cell[c].or2frac[2][2],0.0);

      for (unsigned a                  = 0; a < mdl.atom.size() ; a++)
        if (c                         == mdl.atom[a].crystal)
        {
          vector<double> ortho = mdl.site[mdl.atom[a].nsite].ortho(xtal.cell[c].frac2or);
          if (mdl.atom[a].name.size() == 1)
            fprintf(pdbfile,"HETATM %4d %1s    %2s    %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                a+1,mdl.atom[a].name.substr(0,1).c_str(),
                mdl.atom[a].name.substr(0,2).c_str(),a+1,
                ortho[0],ortho[1],ortho[2],mdl.atom[a].occ,mdl.atom[a].biso);
          else
            fprintf(pdbfile,"HETATM %4d %2s   %2s    %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                a+1,mdl.atom[a].name.substr(0,2).c_str(),
                mdl.atom[a].name.substr(0,2).c_str(),a+1,
                ortho[0],ortho[1],ortho[2],mdl.atom[a].occ,mdl.atom[a].biso);
        }
      fprintf(pdbfile,"END\n");
      fclose(pdbfile);
    }
  }
}

void Likelihood::writexml() const
{
  string xmlname  = output + ".xml";

  FILE *xmlfile   = fopen(xmlname.c_str(),"w");

  fprintf(xmlfile, "<coordinate_data>\n");
  for (unsigned a = 0; a < mdl.atom.size() ; a++)
  {
    fprintf(xmlfile, "\t<atom id=\"%d\">\n",a+1);
    fprintf(xmlfile, "\t\t<atom_id>%s%d</atom_id>\n",
        mdl.atom[a].name.c_str(),a+1);
    fprintf(xmlfile, "\t\t<element_name>%s</element_name>\n",
        mdl.atom[a].name.c_str());
    unsigned s  = mdl.atom[a].nsite;
    fprintf(xmlfile, "\t\t<x>%7.4f</x>\n",
        mdl.site[s].x[0]);
    fprintf(xmlfile, "\t\t<y>%7.4f</y>\n",
        mdl.site[s].x[1]);
    fprintf(xmlfile, "\t\t<z>%7.4f</z>\n",
        mdl.site[s].x[2]);
    fprintf(xmlfile, "\t\t<occupancy>%7.3f</occupancy>\n",
        mdl.atom[a].occ);
    fprintf(xmlfile, "\t\t<bfactor>%7.3f</bfactor>\n",
        mdl.atom[a].biso);
    fprintf(xmlfile, "\t</atom>\n");
  }
  fprintf(xmlfile, "</coordinate_data>\n");
  fprintf(xmlfile, "\n\n");

  fclose(xmlfile);
}

bool Likelihood::hermitianinverse(Matrix &reori, Matrix &imori, Matrix &reinv, 
    Matrix &iminv, double &det, const bool verbose)
{
  // pseudoinverse for a hermitian matrix
  vector<std::complex<double> > cevectors(reori.size());
  evalues.resize(reori.rsize());

  for (unsigned i        = 0; i < reori.size(); i++)
    cevectors[i]         = std::complex<double> (reori.vec(i), imori.vec(i));

  char jobz('V'), uplo('U');
  int n(reori.rsize()), lda(reori.rsize());
  int clwork(n*2 + n*n), clrwork(1 + 5*n + 2*n*n), cliwork(3 + 5*n);
  vector<std::complex<double> > cwork(clwork);
  vector<double> crwork(clrwork, ZERO);
  vector<int> ciwork(cliwork, 0);

  FORTRAN_CALL ( ZHEEVD, zheevd,
      (&jobz, &uplo, &n, &cevectors[0], &lda, &evalues[0], &cwork[0],
       &clwork, &crwork[0], &clrwork,  &ciwork[0], &cliwork, &info),
      (&jobz, &uplo, &n, &cevectors[0], &lda, &evalues[0], &cwork[0],
       &clwork, &crwork[0], &clrwork,  &ciwork[0], &cliwork, &info),
      (&jobz, &uplo, &n, &cevectors[0], &lda, &evalues[0], &cwork[0],
       &clwork, &crwork[0], &clrwork,  &ciwork[0], &cliwork, &info));

  for (unsigned i        = 0; i < reori.rsize(); i++)
    for (unsigned j      = i; j < reori.csize(); j++)
    {
      double recomp(ZERO), recompinv(ZERO);
      double imcomp(ZERO), imcompinv(ZERO);
      for (unsigned k    = 0; k < reori.rsize(); k++)
        if (evalues[k]  >= MINEIG)
        {
          double temp    = (std::real(cevectors[i+reori.rsize()*k])*
              std::real(cevectors[j+reori.rsize()*k]) +
              std::imag(cevectors[i+reori.rsize()*k])*
              std::imag(cevectors[j+reori.rsize()*k]));
          recomp        += temp*evalues[k];
          recompinv     += temp/evalues[k];
          temp           = (std::real(cevectors[i+reori.rsize()*k])*
              std::imag(cevectors[j+reori.rsize()*k]) -
              std::imag(cevectors[i+reori.rsize()*k])*
              std::real(cevectors[j+reori.rsize()*k]));
          imcomp        -= temp*evalues[k];
          imcompinv     -= temp/evalues[k];
        }
      reori(i,j)         = recomp;
      reori(j,i)         = recomp;
      reinv(i,j)         = recompinv;
      reinv(j,i)         = recompinv;
      imori(i,j)         = imcomp;
      imori(j,i)         = -imcomp;
      iminv(i,j)         = imcompinv;
      iminv(j,i)         = -imcompinv;
    }

  det                    = ONE;

  bool filter(false);

  for (unsigned i        = 0; i < reori.rsize(); i++)
  {
    if (evalues[i]       < MINEIG)
      filter             = true;
    det                 *= (evalues[i] < MINEIG) ? MINEIG : evalues[i];
  }

  if (filter && verbose)
  {
    for (unsigned i    = 0; i < reori.rsize(); i++)
      printf("%f ",evalues[i]);
    printf("\n");
  }

  return filter;
}

bool Likelihood::inverse(Matrix &ori, Matrix &inv, double &det)
{
  return inverse_gold(ori, inv, det);
  const int tagMain = Tinverse.start("Main");

  const int n = ori.csize();
  assert(ori.csize() == ori.rsize());
  assert((int)ori.size() == n*n);
  assert(n == 4);

  int errorHandler;
  double lapackWorkspace[16];
  int N = 4;
  int lwork = N*N;
  char chU[] = "U";

  int tag = Tinverse.start("DPOTRF");
  Matrix inv1(ori);
  dpotrf_(chU, &N, inv1.array(), &N, &errorHandler);
  assert(errorHandler >= 0);
  Tinverse.stop(tag);

  if (errorHandler > 0)
  {
    Tinverse.stop(tagMain);
    return inverse_gold(ori, inv, det);
  }

  det = 1.0;
  for (int i = 0; i < n; i++)
    det *= inv1(i,i);
  det *= det;
  assert(det > 0.0);


  tag = Tinverse.start("DPOTRI");
  dpotri_(chU, &N, inv1.array(), &N, &errorHandler);
  assert(0 == errorHandler);
  for (unsigned i        = 0; i < ori.rsize(); i++)
    for (unsigned j      = i; j < ori.csize(); j++)
      inv(j,i) = inv(i,j) = inv1(i,j);
  Tinverse.stop(tag);

  Tinverse.stop(tagMain);
  return false;
}

bool Likelihood::inverse_gold(Matrix &ori, Matrix &inv, double &det)
{
  const int tagMain = TinverseGold.start("Main");

  // pseudoinverse for real symmetric matrix
  evectors.resize(ori.size());
  evalues.resize(ori.rsize());

  for (unsigned i        = 0; i < evectors.size(); i++)
    evectors[i]          = ori.vec(i);

  char jobz('V'), uplo('U');
  int n(ori.rsize()), lda(ori.rsize());

  int tag1 = TinverseGold.start("DSYEVD");
  FORTRAN_CALL ( DSYEVD, dsyevd,
      (&jobz, &uplo, &n, &evectors[0], &lda, &evalues[0], &lawork[0],
       &lwork, &iwork[0], &liwork, &info),
      (&jobz, &uplo, &n, &evectors[0], &lda, &evalues[0], &lawork[0], 
       &lwork, &iwork[0], &liwork, &info),
      (&jobz, &uplo, &n, &evectors[0], &lda, &evalues[0], &lawork[0], 
       &lwork, &iwork[0], &liwork, &info));
  TinverseGold.stop(tag1);

  tag1 = TinverseGold.start("loop1");
  for (unsigned i        = 0; i < ori.rsize(); i++)
    for (unsigned j      = i; j < ori.csize(); j++)
    {
      assert(ori(i,j) == ori(j,i));
      double recomp(ZERO), recompinv(ZERO);
      for (unsigned k    = 0; k < ori.rsize(); k++)
        if (evalues[k]   > MINEIG)
        {
          double temp    = evectors[i+ori.rsize()*k]*evectors[j+ori.rsize()*k];
          recomp        += temp*evalues[k];
          recompinv     += temp/evalues[k];
        }
      ori(i,j)           = recomp;
      ori(j,i)           = recomp;
      inv(i,j)           = recompinv;
      inv(j,i)           = recompinv;
    }
  TinverseGold.stop(tag1);

  det                    = ONE;

  bool filter(false);

  tag1 = TinverseGold.start("loop2");
  for (unsigned i        = 0; i < ori.rsize(); i++)
  {
    if (evalues[i]       < MINEIG)
      filter             = true;
    det                 *= (evalues[i] < MINEIG) ? MINEIG : evalues[i];
  }
  TinverseGold.stop(tag1);

  TinverseGold.stop(tagMain);
  return filter;  
}

void Likelihood::inverse2by2(Matrix &ori, Matrix &inv, double &det)
{
  // ***NSP error when data[0] == data[3] !!
  // returns the positive semi-definite generalized inverse of a 2x2 Hermitian matrix
  // and corrects the original matrix to maintain positive definiteness

  // returns the positive semi-definite generalized inverse of a 2x2 Hermitian matrix
  // and corrects the original matrix to maintain positive definiteness

  double sumsquared  = ori(0,1)*ori(0,1);

  double temp        = sqrt(FOUR*sumsquared + pow(ori(0,0) - ori(1,1), 2));

  // ensure only positive eigenvalues
  double eigenvalue0 = std::max((ori(0,0) + ori(1,1) - temp)/TWO, MINEIG);
  double eigenvalue1 = std::max((ori(0,0) + ori(1,1) + temp)/TWO, MINEIG);

  det                = eigenvalue0*eigenvalue1;

  double mu0         = (ori(0,0) - ori(1,1) - temp)/TWO;
  double mu1         = (ori(0,0) - ori(1,1) + temp)/TWO;

  double arg0(ZERO), arg1(ZERO), invarg0(ZERO), invarg1(ZERO);

  double temp0       = mu0*mu0 + sumsquared;
  if ( (temp0        > EPSILON) && (eigenvalue0 > MINEIG) )
  {
    arg0             = ONE/temp0;
    invarg0          = arg0/eigenvalue0;
    arg0            *= eigenvalue0;
  }

  double temp1       = mu1*mu1 + sumsquared;
  if ( (temp1        > EPSILON) && (eigenvalue1 > MINEIG) )
  {
    arg1             = ONE/temp1;
    invarg1          = arg1/eigenvalue1;
    arg1            *= eigenvalue1;
  }

  inv(0,0)           = mu0*mu0*invarg0 + mu1*mu1*invarg1;
  inv(0,1)           = ori(0,1)*(mu0*invarg0 + mu1*invarg1);
  inv(1,0)           = inv(0,1);

  ori(0,0)           = mu0*mu0*arg0 + mu1*mu1*arg1;
  ori(0,1)          *= mu0*arg0 + mu1*arg1;
  ori(1,0)           = ori(0,1);

  if (sumsquared     < EPSILON)
    if (ori(1,1)     > EPSILON)
      inv(1,1)       = ONE/ori(1,1);
    else
    {
      inv(1,1)       = ZERO;
      ori(1,1)       = ZERO;
    }
  else
  {
    inv(1,1)         = sumsquared*(invarg0 + invarg1);
    ori(1,1)         = sumsquared*(arg0 + arg1);
  }  
}

void Likelihood::hermitianmatrixprod(Matrix &reout, Matrix &imout,
    const Matrix &reain, const Matrix &imain,
    const double *recin, const double *imcin)
{
  // return = A*C*A for hermitian square matrices A,C
  for (unsigned i              = 0; i < reain.rsize(); i++)
    for (unsigned j            = i; j < reain.csize(); j++)
    {
      double recomp(ZERO), imcomp(ZERO);
      for (unsigned k          = 0; k < reain.csize(); k++)
      {
        double recomp1(ZERO), imcomp1(ZERO);
        for (unsigned l        = 0; l < reain.rsize(); l++)
        {
          recomp1             += (reain(i,l)*recin[l+reain.rsize()*k] -
              imain(i,l)*imcin[l+reain.rsize()*k]);
          imcomp1             += (reain(i,l)*imcin[l+reain.rsize()*k] +
              imain(i,l)*recin[l+reain.rsize()*k]);

        }
        recomp                += recomp1*reain(k,j) - imcomp1*imain(k,j);
        imcomp                += recomp1*imain(k,j) + imcomp1*reain(k,j);
      }

      reout(i,j)               = recomp;
      reout(j,i)               = recomp;
      if (i                   == j)
        imout(i,i)             = ZERO;
      else
      {
        imout(i,j)             = imcomp;
        imout(j,i)             = -imcomp;
      }
    }
}

void Likelihood::matrixprod(Matrix &reout, const Matrix &reain, const double *recin)
{
  // return reout = A*C*A for real symmetric square matrices A,C
  for (unsigned i              = 0; i < reain.rsize(); i++)
    for (unsigned j            = i; j < reain.csize(); j++)
    {
      double recomp(ZERO);
      for (unsigned l          = 0; l < reain.rsize(); l++)
        for (unsigned k        = 0; k < reain.csize(); k++)
          recomp              += reain(i,l)*recin[l+reain.rsize()*k]*reain(k,j);
      reout(i,j)               = recomp;
      reout(j,i)               = recomp;
    }
}

void Likelihood::storeinitialcolumns(CMtz::MTZ *MTZIN, vector<float> &fdata,
    const unsigned inc, const unsigned r)
{
  // store initial columns needed for output mtz

  if (allin)
  {
    vector<int> flags(fdata.size());
    float res;
    CMtz::ccp4_lrreff(MTZIN, &res, &fdata[0], &flags[0], (const CMtz::MTZCOL**)
        &colin[0], colin.size(), r+1);
    for (unsigned i           = 0; i < colin.size(); i++)
      if (flags[i])
        fdata[i]              = CCP4::ccp4_nan().f;
  }
  else
  {
    fdata[0]                  = (float) xtal.miller[r][0];
    fdata[1]                  = (float) xtal.miller[r][1];
    fdata[2]                  = (float) xtal.miller[r][2];
  }

  if (xtal.sf[0].use(r)      && (mode != "INTENSITY") && (mode != "DIFFE") && (mode != "GIAC") && (mode != "DELTA") && (mode != "RESCALE") && (mode != "MULT") && (mode != "EVALUES") )
  {
    if ( (target             == "SAD") && (xtal.sf[0].usep(r)) )
    {
      fdata[inc]              = (float) xtal.sf[0].datap[r];
      fdata[inc + 1]          = (float) xtal.sf[0].devp[r];
    }
    else
    {
      fdata[inc]              = (float) xtal.sf[0].datamean(r);
      fdata[inc + 1]          = (float) xtal.sf[0].devmean(r);
    }
  }
}

void Likelihood::setupmtz(CMtz::MTZ *MTZOUT, CMtz::MTZ *MTZIN)
{
  colin.resize(0);
  colout.resize(0);

  CMtz::ccp4_lwtitl(MTZOUT, title.c_str(), 0);
  MTZOUT->refs_in_memory             = 0;

  MTZOUT->mtzsymm.spcgrp             = xtal.sg.number;
  MTZOUT->mtzsymm.nsym               = xtal.sg.NSYM;
  MTZOUT->mtzsymm.nsymp              = xtal.sg.NSYMP;
  MTZOUT->mtzsymm.symtyp             = xtal.sg.lattice;
  CSym::ccp4spg_to_shortname(MTZOUT->mtzsymm.spcgrpname,
      xtal.sg.Gethmname().substr(0,10).c_str());
  strncpy(MTZOUT->mtzsymm.pgname, xtal.sg.Getpointgroup().substr(0,10).c_str(), 11);
  for (unsigned s                    = 0; s < xtal.sg.NSYM; s++ )
  {
    for (unsigned i                  = 0; i < 3; i++ )
      for (unsigned j                = 0; j < 3; j++ )
        MTZOUT->mtzsymm.sym[s][i][j] = xtal.sg.symrot[s][j][i];
    for (unsigned i                  = 0; i < 3; i++ )
      MTZOUT->mtzsymm.sym[s][i][3]   = xtal.sg.symtran[s][i];
  }

  string mtzoutname;
  bool scalepack(false);
  for (unsigned d                    = 0; d < xtal.sf.size(); d++)
    scalepack                        = xtal.sf[d].Getscain().size();

  if ( (getenv("HKLOUT")            != NULL) && !scalepack)
    mtzoutname                       = getenv("HKLOUT");
  else if ( (getenv("HKLOUT2")      != NULL) && scalepack)
    mtzoutname                       = getenv("HKLOUT2");
  else
  {
    mtzoutname                       = "intensities.mtz";
    allin                            = false;
  }

  if (hand                           > 1)
    mtzoutname                      += "-oh";

  MTZOUT->fileout                    = CMtz::MtzOpenForWrite(mtzoutname.c_str());

  if (allin) // ***NSP - this may have to change for the hand problem!!!
    for (int x                       = 0; x < CMtz::MtzNxtal(MTZIN); x++) 
    {
      CMtz::MTZXTAL *xtl             = CMtz::MtzIxtal(MTZIN,x);
      CMtz::MTZXTAL *xo              = CMtz::MtzAddXtal(MTZOUT, xtl->xname, xtl->pname,
          xtl->cell);
      for (int s                     = 0; s < CMtz::MtzNsetsInXtal(xtl); s++) 
      {
        CMtz::MTZSET *set            = CMtz::MtzIsetInXtal(xtl,s);
        CMtz::MTZSET *so             = CMtz::MtzAddDataset(MTZOUT, xo, set->dname,
            set->wavelength);
        for (int c                   = 0; c < CMtz::MtzNcolsInSet(set); c++)
        {
          colin.push_back(CMtz::MtzIcolInSet(set,c));
          colout.push_back(CMtz::MtzAddColumn(MTZOUT, so,
                CMtz::MtzIcolInSet(set,c)->label,
                CMtz::MtzIcolInSet(set,c)->type));
        }
      }
    }

  if (mode                          == "INTENSITY")
  {
    for (unsigned d                  = 0; d < xtal.sf.size(); d++)
    {
      // ***NSP
      string junk                    = "BP3I";

      CMtz::MTZXTAL *xtl             = CMtz::MtzAddXtal(MTZOUT, junk.c_str(),
          junk.c_str(),
          &xtal.Getcell(d)[0]);

      CMtz::MTZSET *set              = CMtz::MtzAddDataset(MTZOUT, xtl,
          xtl->xname,
          0.0);      

      if (!allin &&              (d == 0))
      {
        // h column setup
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "H", "H"));
        // k column setup
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "K", "H"));
        // l column setup
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "L", "H"));
      }

      char datanumber[3];
      sprintf(datanumber,"%1d",d+1);
      string i("I"), sigi("SIGI");

      if (xtal.sf[d].anomalous)
      {
        string label                 = i    + datanumber + "(+)";
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "K"));
        label                        = sigi + datanumber + "(+)";
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "M"));
        label                        = i    + datanumber + "(-)";
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "K"));
        label                        = sigi + datanumber + "(-)";
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "M"));
      }
      else
      {
        string label                 = i    + datanumber;
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "J"));
        label                        = sigi + datanumber;
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "Q"));
      }
    }
  }
  else if (mode                     == "EVALUES")
  {
    // ***NSP - only normalize first data set
    unsigned d(0);

    string junk                      = "BP3I";

    CMtz::MTZXTAL *xtl               = CMtz::MtzAddXtal(MTZOUT, junk.c_str(),
        junk.c_str(),
        &xtal.Getcell(d)[0]);

    CMtz::MTZSET *set                = CMtz::MtzAddDataset(MTZOUT, xtl,
        xtl->xname,
        0.0);      

    if (!allin &&                 (d == 0))
    {
      // h column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "H", "H"));
      // k column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "K", "H"));
      // l column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "L", "H"));
    }

    if (xtal.sf[d].anomalous)
    {
      string label                   = mtzea    +  "(+)";
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "R"));
      label                          = mtzsigea +  "(+)";
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "Q"));
      label                          = mtzea    +  "(-)";
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "R"));
      label                          = mtzsigea +  "(-)";
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, label.c_str(), "Q"));
    }
    else
    {
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzea.c_str(), "R"));
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzsigea.c_str(), "Q"));
    }      
  }
  else if (mode                     == "RESCALE")
  {
    for (unsigned d                  = 0; d < xtal.sf.size(); d++)
    {
      // ***NSP
      string junk                    = "RESCALE";

      CMtz::MTZXTAL *xtl             = CMtz::MtzAddXtal(MTZOUT, junk.c_str(),
          junk.c_str(),
          &xtal.Getcell(d)[0]);

      CMtz::MTZSET *set              = CMtz::MtzAddDataset(MTZOUT, xtl,
          xtl->xname,
          0.0);

      if (!allin &&              (d == 0))
      {
        // h column setup
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "H", "H"));
        // k column setup
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "K", "H"));
        // l column setup
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "L", "H"));
      }

      string scoldatap("G"), scoldevp("L"), scoldata("F"), scoldev("Q");
      if (xtal.sf[d].type           == "INTENSITY")
      {
        scoldatap                    = "K";
        scoldevp                     = "M";
        scoldata                     = "J";
        scoldev                      = "Q";	    
      }

      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, xtal.sf[d].sdata.c_str(), scoldata.c_str()));
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, xtal.sf[d].sdev.c_str(), scoldev.c_str()));

      if (xtal.sf[d].anomalous)
      {
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, xtal.sf[d].sdatap.c_str(), scoldatap.c_str()));
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, xtal.sf[d].sdevp.c_str(), scoldevp.c_str()));
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, xtal.sf[d].sdatam.c_str(), scoldatap.c_str()));
        colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, xtal.sf[d].sdevm.c_str(), scoldevp.c_str()));
      }
    }
  }
  else if ( (mode                   == "DIFFE") || (mode == "GIAC") || (mode == "DELTA") || (mode == "MULT") )
  {
    unsigned d(0);

    // ***NSP
    string junk                      = "AFRO";

    CMtz::MTZXTAL *xtl               = CMtz::MtzAddXtal(MTZOUT, junk.c_str(),
        junk.c_str(),
        &xtal.Getcell(d)[0]);

    CMtz::MTZSET *set                = CMtz::MtzAddDataset(MTZOUT, xtl,
        xtl->xname,
        ZERO);
    if (!allin)
    {
      // h column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "H", "H"));
      // k column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "K", "H"));
      // l column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "L", "H"));
    }

    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzfa.c_str(), "F"));
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzsigfa.c_str(), "Q"));
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzea.c_str(), "R"));
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzsigea.c_str(), "Q"));
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzalpha.c_str(), "P"));
  }
  else
  {
    // ***NSP - I am not sure what I should call the new project/crystal/data set
    // that is to be created.

    string junk                    = "BP3-PHASED";

    if (protocol                  == "PHASECOMB")
      junk                         = "MULTICOMB-PHASED";

    CMtz::MTZXTAL *xtl             = CMtz::MtzAddXtal(MTZOUT,
        junk.c_str(),
        junk.c_str(),
        &xtal.Getcell(0)[0]);

    CMtz::MTZSET *set              = CMtz::MtzAddDataset(MTZOUT, xtl,
        xtl->xname,
        ZERO);
    if (!allin)
    {
      // h column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "H", "H"));
      // k column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "K", "H"));
      // l column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "L", "H"));
    }

    // observed amplitude column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzf.c_str(), "F"));
    // observed sigma column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzsigf.c_str(), "Q"));
    // best amplitude column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzfb.c_str(), "F"));
    // best phase column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzpb.c_str(), "P"));
    // FOM column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzfom.c_str(), "W"));
    // HLA column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzhla.c_str(), "A"));
    // HLB column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzhlb.c_str(), "A"));
    // HLC column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzhlc.c_str(), "A"));
    // HLD column setup
    colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzhld.c_str(), "A"));
    // Difference MAP coefficient
    if ( (target              == "PSAD") || (target == "SAD") ||
        (target              == "MLHL") || (target == "PSD2")  )
    {
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzfdiff.c_str(), "F"));
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzpdiff.c_str(), "P"));
    }

    if ( (target              == "SAD") && (outputhcalc) )
    {
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "HFCALC", "F"));
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "HPCALC", "P"));
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, "HECALC", "R"));
    }

    if (protocol              == "PHASECOMB")
    {
      // best amplitude column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzfcomb.c_str(), "F"));
      // best phase column setup
      colout.push_back(CMtz::MtzAddColumn(MTZOUT, set, mtzpcomb.c_str(), "P"));
    }
  }
}

void Likelihood::outputintensities()
{
  double maxi(ZERO), maxs(ZERO);
  for (unsigned r              = 0; r < xtal.maxref; r++)
    for (unsigned d            = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].use(r)   && (xtal.sf[d].type == "AMPLITUDE"))
      {
        maxi                   = std::max(maxi, (double)xtal.sf[d].datamean(r)*
            xtal.sf[d].datamean(r));
        maxs                   = std::max(maxs, (double)xtal.sf[d].devmean(r)*
            sqrt(FOUR*xtal.sf[d].datamean(r)*
              xtal.sf[d].datamean(r) +
              TWO*xtal.sf[d].devmean(r)*
              xtal.sf[d].devmean(r)));
      }

  double maxis                 = std::max((double)maxi, (double)maxs);
  double scale(ONE);
  if (xtal.sf[0].type         == "AMPLITUDE")
    scale                      = (maxis >= 100000.0) ? 25000.0/maxis : ONE;

  CMtz::MTZ *MTZOUT            = CMtz::MtzMalloc(0,0);
  CMtz::MTZ *MTZIN             = CMtz::MtzGet("HKLIN",0);
  setupmtz(MTZOUT,MTZIN);

  unsigned inc                 = (allin) ? colin.size() : 3;

  unsigned columns(0);
  for (unsigned d              = 0; d < xtal.sf.size(); d++)
    if (xtal.sf[d].anomalous)
      columns                 += 4;
    else
      columns                 += 2;

  for (unsigned r              = 0; r < xtal.maxref; r++)
  {    
    // Default value for columns to be written out is MNF
    vector<float> fdata(columns+inc, CCP4::ccp4_nan().f);

    storeinitialcolumns(MTZIN, fdata, inc, r);

    int c(-1);
    for (unsigned d            = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].anomalous)
      {
        if (xtal.sf[d].usep(r))
        {
          if (xtal.sf[d].type == "AMPLITUDE")
          {
            fdata[inc + ++c]   = scale*xtal.sf[d].datap[r]*xtal.sf[d].datap[r];
            fdata[inc + ++c]   = (scale*xtal.sf[d].devp[r]*
                sqrt(FOUR*xtal.sf[d].datap[r]*xtal.sf[d].datap[r] +
                  TWO*xtal.sf[d].devp[r]*xtal.sf[d].devp[r]));
          }
          else
          {
            fdata[inc + ++c]   = xtal.sf[d].datap[r];
            fdata[inc + ++c]   = xtal.sf[d].devp[r];
          }
        }
        else
          c                   += 2;

        if (xtal.sf[d].usem(r) && !xtal.centric[r])
        {
          if (xtal.sf[d].type == "AMPLITUDE")
          {
            fdata[inc + ++c]   = scale*xtal.sf[d].datam[r]*xtal.sf[d].datam[r];
            fdata[inc + ++c]   = (scale*xtal.sf[d].devm[r]*
                sqrt(FOUR*xtal.sf[d].datam[r]*xtal.sf[d].datam[r] +
                  TWO*xtal.sf[d].devm[r]*xtal.sf[d].devm[r]));
          }
          else
          {
            fdata[inc + ++c]   = xtal.sf[d].datam[r];
            fdata[inc + ++c]   = xtal.sf[d].devm[r];

          }
        }
        else
          c                   += 2;
      }
      else
        if (xtal.sf[d].use(r))
        {
          if (xtal.sf[d].type == "AMPLITUDE")
          {
            fdata[inc + ++c]   = scale*xtal.sf[d].datamean(r)*xtal.sf[d].datamean(r);
            fdata[inc + ++c]   = (scale*xtal.sf[d].devmean(r)*
                sqrt(FOUR*xtal.sf[d].datamean(r)*xtal.sf[d].datamean(r) +
                  TWO*xtal.sf[d].devmean(r)*
                  xtal.sf[d].devmean(r)));
          }
          else
          {
            fdata[inc + ++c]   = xtal.sf[d].datamean(r);
            fdata[inc + ++c]   = xtal.sf[d].devmean(r);	    
          }

        }
        else
          c                   += 2;
    CMtz::ccp4_lwrefl(MTZOUT, &fdata[0], &colout[0], fdata.size(), r+1);
  }

  if (allin)
    CMtz::MtzFree(MTZIN);
  CMtz::MtzPut(MTZOUT, " ");
  CMtz::MtzFree(MTZOUT);
}

void Likelihood::outputhklfile() const
{
  if (xtal.sf.size() > 1)
    Bp3Error("Likelihood::outputhklfile","can only handle one data set at the moment");

  double maxd(ZERO), maxs(ZERO);
  for (unsigned r          = 0; r < xtal.maxref; r++)
    for (unsigned d        = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].use(r))
      {
        maxd               = std::max(maxd, (double)xtal.sf[d].datamean(r));
        maxs               = std::max(maxs, (double)xtal.sf[d].devmean(r));
      }

  double maxds             = std::max((double)maxd, (double)maxs);
  double scale             = (maxds >= 9999.9) ? 9999.9/maxds : ONE;

  FILE *dfile(NULL);
  string filename           = output + ".hkl";
  dfile                     = fopen(filename.c_str(),"w");
  /*
     fprintf(dfile,"TITLE SHELX format from BP3\n");
     fprintf(dfile, "CELL %7.5f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
     1.54178,xtal.cell[0].a, xtal.cell[0].b, xtal.cell[0].c,
     xtal.cell[0].alpha,xtal.cell[0].beta,xtal.cell[0].gamma);
     fprintf(dfile,"ZERR    1.00    0.001   0.001   0.001   0.00   0.00   0.00\n");
     fprintf(dfile,"LATT -7\n");
     fprintf("SYMM -X, Y, -Z\n");
     fprintf(dfile,"HKLF  3\n");
     */
  for (unsigned r          = 0; r < xtal.maxref; r++)
    for (unsigned d        = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].use(r))
        if (xtal.sf[d].anomalous && !xtal.centric[r])
        {
          if (xtal.sf[d].usep(r))
            fprintf(dfile, "%4d%4d%4d%8.2f%8.2f\n",
                xtal.miller[r][0], xtal.miller[r][1], xtal.miller[r][2],
                scale*xtal.sf[d].datap[r], scale*xtal.sf[d].devp[r]);
          if (xtal.sf[d].usem(r))
          {
            int h = (xtal.miller[r][0] == 0) ? 0 : -xtal.miller[r][0];
            int k = (xtal.miller[r][1] == 0) ? 0 : -xtal.miller[r][1];
            int l = (xtal.miller[r][2] == 0) ? 0 : -xtal.miller[r][2];

            fprintf(dfile, "%4d%4d%4d%8.2f%8.2f\n",
                h,k,l,scale*xtal.sf[d].datam[r], scale*xtal.sf[d].devm[r]);
          }
        }
        else
          fprintf(dfile, "%4d%4d%4d%8.2f%8.2f\n",
              xtal.miller[r][0], xtal.miller[r][1], xtal.miller[r][2],
              scale*xtal.sf[d].datap[r], scale*xtal.sf[d].devp[r]);

  fclose(dfile);
}

void Likelihood::outputscafile() const
{
  if (xtal.sf.size() > 1)
    Bp3Error("Likelihood::outputscafile","can only handle one data set at the moment");

  if (xtal.sf[0].type     != "INTENSITY")
    Bp3Error("Likelihood::outputscafile","Intensities must be input");

  double maxd(ZERO), maxs(ZERO);
  for (unsigned r          = 0; r < xtal.maxref; r++)
    for (unsigned d        = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].use(r))
      {
        maxd               = std::max(maxd, (double)xtal.sf[d].datamean(r));
        maxs               = std::max(maxs, (double)xtal.sf[d].devmean(r));
      }

  double maxds             = std::max((double)maxd, (double)maxs);
  double scale             = (maxds >= 99999.9) ? 99999.9/maxds : ONE;

  FILE *dfile(NULL);
  string filename           = output + ".sca";
  dfile                     = fopen(filename.c_str(),"w");
  fprintf(dfile,"    1\n");
  fprintf(dfile," -987\n");

  string name;
  for (unsigned i          = 0; i < xtal.sg.hmname.size(); i++)
    if (!isspace(xtal.sg.hmname[i]))
      name.push_back(tolower(xtal.sg.hmname[i]));

  fprintf(dfile,"%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f %s\n",xtal.cell[0].a,
      xtal.cell[0].b, xtal.cell[0].c, xtal.cell[0].alpha,
      xtal.cell[0].beta,xtal.cell[0].gamma,name.c_str());

  for (unsigned r          = 0; r < xtal.maxref; r++)
    for (unsigned d        = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].use(r))
        if (xtal.sf[d].anomalous && !xtal.centric[r])
        {
          if (xtal.sf[d].anouse(r))
            fprintf(dfile, "%4d%4d%4d%8.1f%8.1f%8.1f%8.1f\n",
                xtal.miller[r][0], xtal.miller[r][1], xtal.miller[r][2],
                scale*xtal.sf[d].datap[r], scale*xtal.sf[d].devp[r],
                scale*xtal.sf[d].datam[r], scale*xtal.sf[d].devm[r]);
          else if (xtal.sf[d].usep(r))
            fprintf(dfile, "%4d%4d%4d%8.1f%8.1f\n",
                xtal.miller[r][0], xtal.miller[r][1], xtal.miller[r][2],
                scale*xtal.sf[d].datap[r], scale*xtal.sf[d].devp[r]);
          else if (xtal.sf[d].usem(r))
            fprintf(dfile, "%4d%4d%4d%8.1f%8.1f\n",
                xtal.miller[r][0], xtal.miller[r][1], xtal.miller[r][2],
                scale*xtal.sf[d].datam[r], scale*xtal.sf[d].devm[r]);
        }
        else
          fprintf(dfile, "%4d%4d%4d%8.1f%8.1f\n",
              xtal.miller[r][0], xtal.miller[r][1], xtal.miller[r][2],
              scale*xtal.sf[d].datap[r], scale*xtal.sf[d].devp[r]);


  fclose(dfile);
}
void Likelihood::writescript(const string program) const
{
  string scriptname         = output + ".sh";

  FILE *sfile               = fopen(scriptname.c_str(),"w");

  fprintf(sfile,"#!/bin/sh\n\n");

  if (command.size())
    fprintf(sfile, "%s HKLIN %s HKLOUT nextrun.mtz << eof > nextrun.log\n\n",command.c_str(), getenv("HKLIN"));
  else
    fprintf(sfile, "%s HKLIN %s HKLOUT nextrun.mtz << eof > nextrun.log\n\n",program.c_str(), getenv("HKLIN"));

  for (unsigned s           = 0; s < mdl.site.size() ; s++) 	
  {
    fprintf(sfile,"SITE %3d %5.3f %5.3f %5.3f  ", s + 1,
        mdl.site[s].x[0], mdl.site[s].x[1], mdl.site[s].x[2]);

    if (!mdl.site[s].Getrefinex(0) || !mdl.site[s].Getrefinex(1) ||
        !mdl.site[s].Getrefinex(2))
    {
      fprintf(sfile, " NOREf ");
      if (!mdl.site[s].Getrefinex(0))
        fprintf(sfile, " X ");
      if (!mdl.site[s].Getrefinex(1))
        fprintf(sfile, " Y ");
      if (!mdl.site[s].Getrefinex(2))
        fprintf(sfile, " Z ");
    }
    fprintf(sfile, "\n");
  }

  fprintf(sfile, "\n\n");

  for (unsigned c           = 0; c < xtal.cell.size(); c++)
  {
    fprintf(sfile, "XTAL %s\n", xtal.cell[c].name.c_str());
    int w                   = -1;
    fprintf(sfile, "  CELL %8.3f %8.3f %8.3f %7.3f %7.3f %7.3f\n",
        xtal.cell[c].a, xtal.cell[c].b, xtal.cell[c].c,
        xtal.cell[c].alpha, xtal.cell[c].beta,
        xtal.cell[c].gamma);

    for (unsigned a         = 0; a < mdl.atom.size() ; a++)
      if (c                == mdl.atom[a].crystal)
      {
        fprintf(sfile, "  ATOM %s SITE %u\n", mdl.atom[a].name.c_str(),
            mdl.atom[a].nsite + 1);
        fprintf(sfile, "    OCCU %7.3f", mdl.atom[a].occ);
        if (!mdl.atom[a].Getrefineocc())
          fprintf(sfile, " NOREf");
        fprintf(sfile, "\n");
        if (mdl.atom[a].isotropic)
          fprintf(sfile, "    BISO %7.3f", mdl.atom[a].biso);
        else
        {
          fprintf(sfile, "    UANO ");
          for (unsigned b   = 0; b < 6; b++)
            fprintf(sfile, "%7.3f ",mdl.atom[a].uaniso[b]);
        }
        if (!mdl.atom[a].Getrefinebfac())
          fprintf(sfile, " NOREf");
        fprintf(sfile, "\n");
      }

    for (unsigned d         = 0; d < xtal.sf.size()   ; d++)
      if (c                == xtal.sf[d].nxtal)
      {
        w++;
        fprintf(sfile, "  DNAMe %s\n", xtal.sf[d].name.c_str());
        if (xtal.sf[d].anomalous)
          fprintf(sfile, "    COLUmn  F+= %s  SF+= %s  F-= %s  SF-= %s\n",
              xtal.sf[d].sdatap.c_str(), xtal.sf[d].sdevp.c_str(),
              xtal.sf[d].sdatam.c_str(), xtal.sf[d].sdevm.c_str());
        else
          fprintf(sfile, "    COLUmn  F= %s  SF= %s\n",
              xtal.sf[d].sdatap.c_str(),  xtal.sf[d].sdevp.c_str());

        fprintf(sfile, "    RESOlution %7.2f %7.2f\n", 
            xtal.sf[d].hires, xtal.sf[d].lowres);
        fprintf(sfile, "    BINS %u\n", xtal.sf[d].nbins);
        if ( (d            != 0) && (target == "UNCO") )
        {
          fprintf(sfile, "    ISOE ");
          for (unsigned s   = 0; s < xtal.sf[d].nbins; s++)
            fprintf(sfile, "%9.6f ",xtal.sf[d].dluz[s]);
          fprintf(sfile, "\n");
        }
        if (xtal.sf[d].refinee)
        {
          fprintf(sfile, "    CORE ");
          for (unsigned s   = 0; s < xtal.sf[d].nbins; s++)
            fprintf(sfile, "%9.6f ",xtal.sf[d].eluz[s]);
          fprintf(sfile, "\n");
        }

        if (xtal.sf[d].anomalous)
        {
          if (target       == "UNCO")
          {
            fprintf(sfile, "    ANOE ");
            for (unsigned s = 0; s < xtal.sf[d].nbins; s++)
              fprintf(sfile, "%9.6f ",xtal.sf[d].adluz[s]);
            fprintf(sfile, "\n");
          }
          else if (target  == "SAD")
          {
            fprintf(sfile, "    SDLU ");
            for (unsigned s = 0; s < xtal.sf[d].nbins; s++)
              fprintf(sfile, "%9.6f ",xtal.sf[d].sdluz[s]);
            fprintf(sfile, "\n");
          }
        }
        for (unsigned f     = 0; f < mdl.form.size(); f++)
          if (d            == xtal.dataset[mdl.form[f].Getcrystal()][w])
            fprintf(sfile, "    FORM  %s FP= %7.3f  FPP= %7.3f\n",
                mdl.form[f].Getname().c_str(), mdl.form[f].Getfp(w), 
                mdl.form[f].Getfpp(w));
      }
    fprintf(sfile, "\n\n");
  }
  fprintf(sfile,"TARGet %s\n",target.c_str());
  if ( (target             == "MSRS") || (target == "PSAD") || (target == "PSD2") )
  {
    unsigned d(0);
    if (target == "MSRS")
      d++;

    for (unsigned p         = 0; p < xtal.pluz.size(); p++)
    {
      fprintf(sfile,"PLUZzati %d  ",p+1);
      for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
        fprintf(sfile, "%9.6f ",xtal.pluz[p][s]);
      if (!xtal.refinep[p])
        fprintf(sfile, "NOREf ");
      fprintf(sfile, "\n");
    }
  }
  fprintf(sfile, "\n");
  fprintf(sfile, "# Refine all parameters\n");
  fprintf(sfile, "REFAll\n");
  fprintf(sfile, "\n\n");
  fprintf(sfile, "eof\n");  
  fclose(sfile);
}

void Likelihood::printfom() const
{
  printf(" $TABLE: Figure of merit vs. Resolution-centric, acentric, overall:\n");
  printf("$GRAPHS: Figure of merit vs. Res :N:4,7,9,11:\n");
  printf("       : Figure of merit vs. STOL2 :N:5,7,9,11:\n");

  double totalcfom(ZERO);
  unsigned totalcnshl(0);

  double totalafom(ZERO);
  unsigned totalanshl(0);

  printf("$$\n Bin   LoRes  HiRes    Res   STOL2    Refls  Centr    Refls  Acentr   Refls   All  $$\n");

  printf("$$\n");

  for (unsigned s = 0; s < cfom.size(); s++)
  {
    double cfbin  = (cnshl[0][s]) ? cfom[s]/(double)(cnshl[0][s]) : ZERO;
    double afbin  = (anshl[0][s]) ? afom[s]/(double)(anshl[0][s]) : ZERO;
    double tfbin  = ( (anshl[0][s]+cnshl[0][s]) ?
        (afom[s] + cfom[s])/
        (double)(anshl[0][s] + cnshl[0][s]) : ZERO);

    printf(" %2d  %6.2f  %6.2f %6.2f  %7.5f  %5d  %7.5f  %5d  %7.5f  %5d  %7.5f\n",
        s+1, xtal.sf[0].Getlores(s), xtal.sf[0].Gethires(s), xtal.sf[0].Getavres(s),xtal.sf[0].astolsq[s],
        cnshl[0][s], cfbin, anshl[0][s], afbin, anshl[0][s] + cnshl[0][s], tfbin);

    totalcfom    += cfom[s];
    totalcnshl   += cnshl[0][s];
    totalafom    += afom[s];
    totalanshl   += anshl[0][s];
  }
  double totcfom  = (totalcnshl) ? totalcfom/(double)(totalcnshl) : ZERO;
  double totafom  = (totalanshl) ? totalafom/(double)(totalanshl) : ZERO;
  double totfom   = ( (totalanshl+ totalcnshl) ?
      (totalafom+totalcfom)/(double)(totalanshl+ totalcnshl) : ZERO);

  printf("$$\nTOTAL\n");
  printf("     %6.2f  %6.2f                  %5d  %7.5f  %5d  %7.5f  %5d  %7.5f\n\n",
      xtal.sf[0].Getlores(0), xtal.sf[0].Gethires(xtal.sf[0].Getnbins()-1),
      totalcnshl, totcfom, totalanshl, totafom, totalanshl+totalcnshl, totfom);

  // printf("\n\nThe overall FOM is %.3f\n\n",totfom);
  fflush(stdout);
}
