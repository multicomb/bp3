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

#include <math.h>
#include "gauss.h"
#include "bp3likelihood.h"
#include "matrix.h"
#include "csymlib.h"

#include "mytimer.h"
extern Timer Tgrad1dim, Tgrad;

  Bp3likelihood::Bp3likelihood(Model &model, Crystal &crystal)
: Likelihood(model, crystal)
{
  // defaults
  likelihood      = ZERO;
  mode            = "REFINE";
  //  target          = "UNCO";
  protocol        = "PHASING";
  title           = "Phasing from BP3";
  threshold       = FOUR;
  difftol         = 1.0e-7;
  verbose         = cd             = 0;
  shelthres       = 50;
  refineocc       = true;       
  refinescale     = refinexyz      = refinebfac  = refineluzzati = false;
  refinefp        = refinefpp      = refinenumb  = false;
  minvar          = EPSILON;
  allin           = stats          = sheldrick   = refineall     = false;
  interpolate     = false;
  output          = "heavy";
  // default size for Gaussian integration parameters
  cweight.resize(3);
  pweight.resize(30);
  afweight.resize(3);
  sadweight.resize(35);
}


double Bp3likelihood::gradient1dim (double &dfunc, vector<double> &grd,
    const double stepsize,
    const vector<double> &pars,
    const vector<double> &direction)
{ 
  std::string tag;

  tag = Tgrad1dim.start("shift");
  shift(stepsize,pars,direction); 
  Tgrad1dim.stop(tag);

  tag = Tgrad1dim.start("directsfcalc");
  directsfcalc(false,updatesigmah); 
  Tgrad1dim.stop(tag);

  tag = Tgrad1dim.start("grad");
  double func(grad());
  Tgrad1dim.stop(tag);

  tag = Tgrad1dim.start("Getgradient");
  grd = Getgradient();
  Tgrad1dim.stop(tag);

  tag = Tgrad1dim.start("dotprod");
  dfunc = dot_prod(grd,direction); 
  Tgrad1dim.stop(tag);
  return func;
}
double Bp3likelihood::grad(const bool check, const bool outputmtz)
{
  std::string tag;
  // Evaluates the function and the gradient

  // Calculation derivatives of requested function
  double val(ZERO);

  if (target      == "UNCO")
  {
    tag = Tgrad.start("UNCO");
    val            = uncorrelatedgradient(check,outputmtz);
    Tgrad.stop(tag);
  }
  else if (target == "SAD")
  {
    tag = Tgrad.start("SAD");
    val            = sadgradient(check,outputmtz);
    Tgrad.stop(tag);
  }
  else if (target == "MAD")
  {
    tag = Tgrad.start("MAD");
    val            = uncorrelatedgradient(check,outputmtz);
    Tgrad.stop(tag);
  }
  else if (target == "MSRS")
  {
    tag = Tgrad.start("MSRS");
    val            = multsirasgradient(check,outputmtz);
    Tgrad.stop(tag);
  }
  else if (target == "PSAD")
  {
    tag = Tgrad.start("PSAD");
    val            = pavolsadgradient(check,outputmtz);
    Tgrad.stop(tag);
  }

  // The above routine only calculate d{Target}/d{Structure Factor}.
  // The routine below calculates d{Target}/d{atom parameter}, if we
  // are refining atomic parameters.
  //
  tag = Tgrad.start("directsfcalc");
  if (refineocc || refinexyz || refinebfac || (mode == "DIFF"))
    directsfcalc(true,updatesigmah);
  Tgrad.stop(tag);

  return val;
}

double Bp3likelihood::function(const bool check)
{
  double val(ZERO);

  if (target      == "UNCO")
    val            = uncorrelatedgradient(check);
  else if (target == "SAD")
    val            = sadgradient(check);    
  else if (target == "MAD")
    val            = uncorrelatedgradient(check);
  else if (target == "MSRS")
    val            = multsirasfunction(check);
  else if (target == "PSAD")
    val            = pavolsadfunction(check);

  if (refineocc || refinexyz || refinebfac || (mode == "DIFF"))
    directsfcalc(false,updatesigmah);

  return val;
}

void Bp3likelihood::setup(const bool first, const bool refine, const bool outputmtz)
{
  // Determine how many isomorphism and atomic parameters there are  
  // to refine and allocate memory for variables storing gradients.

  if (mode                     == "CHECK")
  {
    if (xtal.sf[0].anomalous)
    {
      target                    = "SAD";
      cd                        = 0;

      if (xtal.sf.size()        > 1)
      {
        vector<double> totaldano(xtal.sf.size(), ZERO);
        vector<double> total(xtal.sf.size(), ZERO);	
        for (unsigned d         = 0; d < xtal.sf.size(); d++)
          for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
          {
            totaldano[d]      += xtal.sf[d].danoovers[s]*xtal.sf[d].anonshl[s];
            total[d]          += (double)xtal.sf[d].anonshl[s];
          }

        for (unsigned d         = 0; d < xtal.sf.size(); d++)
          if (totaldano[d]/total[d] > totaldano[cd]/total[cd])
            cd                  = d;
      }

      for (unsigned d           = 0; d < xtal.sf.size(); d++)
        if (d                  != cd)
          xtal.sf[d].refinesd   = false;
    }    
  }
  else if (xtal.mad)  // ***NSP
  {
    target                      = "UNCO";
    if (!refine)
      madphasingdefaults();
  }

  if (xtal.sf.size()           == 1)
  {
    if (xtal.sf[0].anomalous && !target.size())
      target                    = "SAD";
  }
  else if (!target.size())
    target                      = "UNCO";

  resize(outputmtz);

  if (target                   == "MSRS")
    if (refinexyz && refineocc && refinebfac)
      xtal.refinep[2]           = xtal.refinep[3] = true;
    else
      xtal.refinep[2]           = xtal.refinep[3] = false;

  // Setup Gaussian variables
  Setgauss();

  if (!xtal.mad && (target     == "UNCO") )
    for (unsigned d             = 1; d < xtal.sf.size(); d++)
      for (unsigned s           = 0; s < xtal.sf[d].nbins; s++)
        if (xtal.sf[d].var(ONE, s) <= ONE && (xtal.sf[d].sigman[s] > xtal.sf[d].sigmah[s]))
          xtal.sf[d].dluz[s]    = 0.99*sqrt((xtal.sf[d].sigman[s] -
                xtal.sf[d].sigmah[s])/
              xtal.sf[d].sigmanref[s]);

  if (first && (target         == "MSRS") && !xtal.userpluz[1])
    for (unsigned s             = 0; s < xtal.sf[1].nbins; s++)
      xtal.pluz[1][s]           = xtal.sf[1].dluz[s];


  if (first && (target         == "PSAD") && !xtal.userpluz[0])
    for (unsigned s             = 0; s < xtal.sf[0].nbins; s++)
      xtal.pluz[0][s]           = 0.5;

  if (refine)
  {
    unsigned nscalepars(0);
    unsigned nluzzatipars(0);
    unsigned natmpars(0);

    for (unsigned d             = 0; d < xtal.sf.size(); d++)
    {
      if (target               == "UNCO")
      {
        if (xtal.sf[d].refinek && refinescale)
          nscalepars++;
        if (xtal.sf[d].refineb && refinescale)
          nscalepars++;

        if (xtal.sf[d].refined && refineluzzati)
          nluzzatipars         += xtal.sf[d].nbins;
        if (xtal.sf[d].refinee && refineluzzati)
          nluzzatipars         += xtal.sf[d].nbins;             
        if (xtal.sf[d].refinead && refineluzzati)
          nluzzatipars         += xtal.sf[d].nbins;
      }
      else if ( (target        == "SAD") || (target == "MAD") )
        if (xtal.sf[d].refinesd && refineluzzati)
          nluzzatipars         += xtal.sf[d].nbins;      
    }

    if ( target                == "MSRS")
      for (unsigned p           = 0; p < xtal.pluz.size(); p++)
        if (xtal.refinep[p]    && refineluzzati)
          nluzzatipars         += xtal.sf[1].nbins;

    if ( target                == "PSAD")
      for (unsigned p           = 0; p < xtal.pluz.size(); p++)
        if (xtal.refinep[p]    && refineluzzati)
          nluzzatipars         += xtal.sf[0].nbins;

    // Atom parameters
    if (refinexyz)
      for (unsigned s           = 0; s < mdl.site.size(); s++)
        for (unsigned k         = 0; k < 3              ; k++)
          if (mdl.site[s].Getrefinex(k))
            natmpars++;

    for (unsigned a             = 0; a < mdl.atom.size()  ; a++)
    {
      if (mdl.atom[a].Getrefinebfac() && refinebfac)
        if (mdl.atom[a].Getisotropic() )
          natmpars++;
        else  
          natmpars             += 6;
      if (mdl.atom[a].Getrefineocc() && refineocc)
        natmpars++;      
    }

    for (unsigned f             = 0; f < mdl.form.size()  ; f++)
      for (unsigned w           = 0; w < mdl.form[f].Getnwave(); w++)
      {
        if (mdl.form[f].Getrefinefp(w) && refinefp)
          natmpars++;
        if (mdl.form[f].Getrefinefpp(w) && refinefpp)
          natmpars++;
      }

    npars                       = natmpars + nluzzatipars + nscalepars;
    printf("Number of parameters being refined: %u\n", npars);
    printf("Consisting of...\n");
    if (nscalepars)
      printf("Scale parameters: %u\n", nscalepars);
    if (nluzzatipars)
      printf("Luzzati parameters: %u\n", nluzzatipars);
    if (natmpars)
      printf("Atomic parameters: %u\n\n", natmpars);

    Settype();

    // ensure positive definiteness
    if ( (target               == "MSRS") || (target == "PSAD") )
    {
      double oldval(ZERO);
      for (unsigned i           = 0; i < 10; i++)
      {
        double val(ZERO);
        if (target             == "MSRS")
          val                   = multsirasgradient();
        else if (target        == "PSAD")
          val                   = pavolsadgradient();
        if (fabs(val - oldval)  < DSMALL)
          break;
        oldval                  = val;
      }
    }

    // check for reflections with low-likelihood
    grad(true,outputmtz);
  }

  /*
     if ((xtal.sf[0].sigmah[0]    <= ZERO) && xtal.mad)
     Bp3Error("Bp3likelihood::setup", "define the native data set last in mad");
     */

  // xtal.print();
}

unsigned Bp3likelihood::Setiterations() const
{

  unsigned iterations(6);

  if (mode               != "CHECK" )
  { 
    unsigned wavelengths(0);

    for (unsigned d     = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].nwave)
        wavelengths       = std::max(xtal.sf[d].nwave+1, wavelengths);

    if (wavelengths      == 2)
      iterations          = 10;
    else if (wavelengths == 3)
      iterations          = 15;
    else if (wavelengths == 4)
      iterations          = 20;

  }

  printf("Setting the number of iterations to %u\n\n", iterations);

  return iterations;
}

void Bp3likelihood::checkatoms()
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

void Bp3likelihood::madphasingdefaults()
{
  if (xtal.Getmad())
    xtal.Setonlyacentrics(false).Setonlycentrics(false);
}

Bp3likelihood &Bp3likelihood::Setgauss()
{
  // Setup Gaussian variables

  if ( (target          == "UNCO") && (!gsin.size()))
  {
    Gauss gausscentrichermite("Hermite", cweight.size());
    cnode                = gausscentrichermite.Getnode();
    cweight              = gausscentrichermite.Getweight();

    for (unsigned i      = 0; i < cweight.size(); i++)
      cweight[i]        *= exp(cnode[i]*cnode[i]);

    if (verbose          > 1)
      gausscentrichermite.print();

    Gauss gausslegendre("Legendre", pweight.size());
    vector<double> pnode = gausslegendre.Getnode();
    pweight              = gausslegendre.Getweight();

    if (verbose          > 1)
      gausslegendre.print();

    gsin.resize(pweight.size());
    gcos.resize(pweight.size());  
    // Precompute sin, cos                                             
    for (unsigned i      = 0; i < pweight.size(); i++)
    {
      gcos[i]            = tab.Cos(PI*pnode[i]);
      gsin[i]            = tab.Sin_charged(PI*pnode[i]);
    }

    Gauss gausshermite("Hermite", afweight.size());
    afnode               = gausshermite.Getnode();
    afweight             = gausshermite.Getweight();

    for (unsigned i      = 0; i < afweight.size(); i++)
      afweight[i]       *= exp(afnode[i]*afnode[i]);

    if (verbose          > 1)
      gausshermite.print();

    sadweight.resize(0);
  }
  else if ( ( (target    == "SAD") || (target == "MAD")) && !sadsin.size())
  {
    // Sad function
    Gauss gausslegendre("Legendre", sadweight.size());
    vector<double> snode = gausslegendre.Getnode();
    sadweight            = gausslegendre.Getweight();

    sadsin.resize(sadweight.size());
    sadcos.resize(sadweight.size());  
    // Precompute sin, cos                                             
    for (unsigned i      = 0; i < sadweight.size(); i++)
    {
      sadcos[i]          = tab.Cos(PI*snode[i]);
      sadsin[i]          = tab.Sin_charged(PI*snode[i]);
    }
    cweight.resize(0);
    pweight.resize(0);
    afweight.resize(0);

    if (verbose          > 1)
      gausslegendre.print();
  }
  return *this;
}

void Bp3likelihood::otherhand()
{
  unsigned enantinumber(xtal.sg.enantiomorph());
  output   += "-oh";
  hand      = 2;
  if (enantinumber)
    xtal.sg = Spacegroup(enantinumber);
  mdl.inverthand(enantinumber);
}

void Bp3likelihood::finitedifftest()
{
  // Perform finite difference tests on various derivatives calculated:
  double likelihood2;

  if (target               == "UNCO")
  {
    printf("\nSCALE PARAMETERS\n");

    double value            = grad();

    for (unsigned d         = 0; d < xtal.sf.size(); d++)
    {
      printf("dLdkscale   = %f\n", dLdkscale[d]);
      printf("dLdbscale   = %f\n", dLdbscale[d]);
    }

    unsigned de             = xtal.sf.size()-1;

    xtal.sf[de].biso       += difftol;

    likelihood2             = function();

    printf("Likelihood from function = %f\n", likelihood2);
    printf("Likelihood from gradient = %f\n", value);
    printf("Finite diff              = %f\n", (likelihood2 - value)/difftol);

    value                   = grad();

    for (unsigned d         = 0; d < xtal.sf.size(); d++)
    {
      printf("dLdkscale   = %f\n", dLdkscale[d]);
      printf("dLdbscale   = %f\n", dLdbscale[d]);
    }

    xtal.sf[de].kscale     += difftol;

    likelihood2 = function();    

    printf("Likelihood from function = %f\n", likelihood2);
    printf("Likelihood from gradient = %f\n", value);
    printf("Finite diff              = %f\n", (likelihood2 - value)/difftol);

    printf("\nNON-ISOMORPHISM PARAMETERS\n");    

    value                   = grad();

    for (unsigned d         = 0; d < xtal.sf.size()    ; d++)
    {
      for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
        printf("GradDluz[%u][%u]= %f \n",d,s,dLddluz[d][s]);
      printf("\n");
      for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
        printf("GradEluz[%u][%u]= %f \n",d,s,dLdeluz[d][s]);
      printf("\n");
    }

    unsigned sh             = xtal.sf[de].nbins -1;

    if (xtal.sf[de].refined)
      xtal.sf[de].dluz[sh] += difftol;
    else if (xtal.sf[de].refinee)
      xtal.sf[de].eluz[sh] += difftol;

    likelihood2             = function();

    printf("Likelihood from function = %f\n", likelihood2);
    printf("Likelihood from gradient = %f\n", value);
    printf("Finite diff              = %f\n", (likelihood2 - value)/difftol);

    if (xtal.sf[de].anomalous)
    {
      value                 = grad();

      for (unsigned d       = 0; d < xtal.sf.size()  ; d++)
      {
        for (unsigned s     = 0; s < xtal.sf[d].nbins; s++)
          printf("GradADluz[%u][%u]= %f \n",d,s,dLdadluz[d][s]);
        printf("\n");
      }

      xtal.sf[de].adluz[0] += difftol;

      likelihood2           = function();

      printf("Likelihood from function = %f\n", likelihood2);
      printf("Likelihood from gradient = %f\n", value);
      printf("Finite diff              = %f\n", (likelihood2 - value)/difftol);
    }
  }
  else if ( (target        == "SAD") || (target == "MAD") )
  {
    double value            = grad();

    for (unsigned d         = 0; d < xtal.sf.size()    ; d++)
    {
      for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
        printf("GradSDluz[%u][%u]= %f \n",d,s,dLdsdluz[d][s]);
      printf("\n");
    }

    unsigned s              =xtal.sf[0].nbins -1;

    if (xtal.sf[0].refinesd)
      xtal.sf[0].sdluz[s]  -= difftol;

    likelihood2             = function();

    printf("Likelihood from function = %f\n", likelihood2);
    printf("Likelihood from gradient = %f\n", value);
    printf("Finite diff              = %f\n", (value - likelihood2)/difftol);
  }
  else if ( target         == "MSRS" || target         == "PSAD" )
  {
    unsigned d(0);
    if (target             == "MSRS")
      d++;

    double value            = grad();

    for (unsigned p         = 0; p < xtal.pluz.size()    ; p++)
    {
      for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
        printf("GradPluz[%u][%u]= %f \n",p,s,dLdpluz[p][s]);
      printf("\n");
    }

    unsigned s              = 0;

    xtal.pluz[1][s]        -= difftol;

    likelihood2             = function();

    printf("Likelihood from function = %f\n", likelihood2);
    printf("Likelihood from gradient = %f\n", value);
    printf("Finite diff (PLuz[1][0]) = %f\n", (value - likelihood2)/difftol);

    value                   = grad();

    for (unsigned p         = 0; p < xtal.pluz.size()    ; p++)
    {
      for (unsigned s       = 0; s < xtal.sf[d].nbins; s++)
        printf("GradPluz[%u][%u]= %f \n",p,s,dLdpluz[p][s]);
      printf("\n");
    }

    s                       = 0;

    xtal.pluz[0][s]        -= difftol;

    likelihood2             = function();

    printf("Likelihood from function = %f\n", likelihood2);
    printf("Likelihood from gradient = %f\n", value);
    printf("Finite diff (PLuz[0][0]) = %f\n", (value - likelihood2)/difftol);

  }

  directsfcalc(false,updatesigmah);  

  printf("\nATOMIC PARAMETERS\n");    

  unsigned id1            = mdl.site.size()-1;

  double value            = grad();

  for (unsigned i         = 0; i < mdl.site.size(); i++)
  {
    printf(" Grad[%u][X]   = %f\n", i, dLdx[i][0]);
    printf(" Grad[%u][Y]   = %f\n", i, dLdx[i][1]);
    printf(" Grad[%u][Z]   = %f\n", i, dLdx[i][2]);
  }

  mdl.site[id1].Setx(2,mdl.site[id1].x[2] + difftol);

  directsfcalc(false,updatesigmah);

  likelihood2             = function();

  printf("Likelihood from function = %f\n", likelihood2);
  printf("Likelihood from gradient = %f\n", value);
  printf("Finite diff              = %f\n", (likelihood2 - value)/difftol);

  value                   = grad();

  for (unsigned i         = 0; i < mdl.site.size(); i++)
  {
    printf(" Grad[%u][X]   = %f\n", i, dLdx[i][0]);
    printf(" Grad[%u][Y]   = %f\n", i, dLdx[i][1]);
    printf(" Grad[%u][Z]   = %f\n", i, dLdx[i][2]);
  }

  mdl.site[id1].Setx(0,mdl.site[id1].x[0] + difftol);

  directsfcalc(false,updatesigmah);

  likelihood2             = function();

  printf("Likelihood from function = %f\n", likelihood2);
  printf("Likelihood from gradient = %f\n", value);
  printf("Finite diff              = %f\n", (likelihood2 - value)/difftol); 

  value                   = grad();

  for (unsigned i         = 0; i < mdl.atom.size(); i++)
  {
    printf(" Grad[%u][O]   = %f\n",i,dLdocc[i]);
    if (mdl.atom[id1].isotropic)
      printf(" Grad[%u][B]   = %f\n",i, dLdbiso[i]);
    else
    {
      printf(" Grad[%u][U1]   = %f\n",i, dLduaniso[i][0]);
      printf(" Grad[%u][U2]   = %f\n",i, dLduaniso[i][1]);
      printf(" Grad[%u][U3]   = %f\n",i, dLduaniso[i][2]);
      printf(" Grad[%u][U4]   = %f\n",i, dLduaniso[i][3]);
      printf(" Grad[%u][U5]   = %f\n",i, dLduaniso[i][4]);
      printf(" Grad[%u][U6]   = %f\n",i, dLduaniso[i][5]);
    }
    /*
       printf(" Grad[" << i << "][f']  = " << dLdfp[i]   << "n";
       printf(" Grad[" << i << "][f''] = " << dLdfpp[i]  << "n";
       */
  }

  if (mdl.atom[id1].isotropic)
    mdl.atom[id1].Setbiso(mdl.atom[id1].biso + difftol);
  else
    mdl.atom[id1].Setuaniso(mdl.atom[id1].uaniso[0] + difftol, 0);

  directsfcalc(false,updatesigmah);

  likelihood2             = function();

  printf("Likelihood from function = %f\n", likelihood2);
  printf("Likelihood from gradient = %f\n", value);
  printf("Finite diff              = %f\n", (likelihood2 - value)/difftol);

  value                   = grad();

  for (unsigned i         = 0; i < mdl.atom.size(); i++)
    printf(" Grad[%u][O] = %f\n",i,dLdocc[i]);

  mdl.atom[id1].Setocc(mdl.atom[id1].occ+difftol);

  directsfcalc(false,updatesigmah);

  likelihood2             = function();

  printf("Likelihood from function = %f\n", likelihood2);
  printf("Likelihood from gradient = %f\n", value);
  printf("Finite diff              = %f\n", (likelihood2 - value)/difftol);  
}

void Bp3likelihood::initstat()
{
  aisolofc.resize(xtal.sf.size());
  cisolofc.resize(xtal.sf.size());
  anolofc.resize(xtal.sf.size());
  asumfph.resize(xtal.sf.size());
  csumfph.resize(xtal.sf.size());
  anosumfph.resize(xtal.sf.size());
  asumfh.resize(xtal.sf.size());
  csumfh.resize(xtal.sf.size());
  sumimfh.resize(xtal.sf.size());
  asumdiff.resize(xtal.sf.size());
  csumdiff.resize(xtal.sf.size());
  sanonshl.resize(xtal.sf.size());
  sumdano.resize(xtal.sf.size());

  for (unsigned d   = 0; d < xtal.sf.size(); d++)
  {
    aisolofc[d].resize(xtal.sf[d].nbins, ZERO);
    cisolofc[d].resize(xtal.sf[d].nbins, ZERO);
    anolofc[d].resize(xtal.sf[d].nbins, ZERO);
    asumfph[d].resize(xtal.sf[d].nbins, ZERO);
    csumfph[d].resize(xtal.sf[d].nbins, ZERO);
    anosumfph[d].resize(xtal.sf[d].nbins, ZERO);
    asumfh[d].resize(xtal.sf[d].nbins, ZERO);
    csumfh[d].resize(xtal.sf[d].nbins, ZERO);
    sumimfh[d].resize(xtal.sf[d].nbins, ZERO);
    asumdiff[d].resize(xtal.sf[d].nbins, ZERO);
    csumdiff[d].resize(xtal.sf[d].nbins, ZERO);
    sanonshl[d].resize(xtal.sf[d].nbins, 0);
    sumdano[d].resize(xtal.sf[d].nbins, ZERO);
  }
}

void Bp3likelihood::printstat() const
{
  for (unsigned d          = 0; d < xtal.sf.size(); d++)
  {
    vector<double> aphasepower(xtal.sf[d].nbins, ZERO);
    vector<double> cphasepower(xtal.sf[d].nbins, ZERO);
    vector<double> arkraut(xtal.sf[d].nbins, ZERO);
    vector<double> crkraut(xtal.sf[d].nbins, ZERO);
    vector<double> arcullis(xtal.sf[d].nbins, ZERO);
    vector<double> crcullis(xtal.sf[d].nbins, ZERO);
    vector<double> anophasepower(xtal.sf[d].nbins, ZERO);
    vector<double> anorkraut(xtal.sf[d].nbins, ZERO);
    vector<double> anorcullis(xtal.sf[d].nbins, ZERO);

    for (unsigned s        = 0; s < xtal.sf[d].nbins; s++)
    {
      if (aisolofc[d][s])      
        aphasepower[s]     = asumfh[d][s]/aisolofc[d][s];
      if (cisolofc[d][s])
        cphasepower[s]     = csumfh[d][s]/cisolofc[d][s];
      if (asumfph[d][s])
        arkraut[s]         = aisolofc[d][s]/asumfph[d][s];
      if (csumfph[d][s])
        crkraut[s]         = cisolofc[d][s]/csumfph[d][s];
      if (asumdiff[d][s])
        arcullis[s]        = aisolofc[d][s]/asumdiff[d][s];
      if (csumdiff[d][s])
        crcullis[s]        = cisolofc[d][s]/csumdiff[d][s];
      if (xtal.sf[d].anomalous)
      {
        if (anolofc[d][s])
          anophasepower[s] = sumimfh[d][s]/anolofc[d][s];
        if (anosumfph[d][s])
          anorkraut[s]     = anolofc[d][s]/anosumfph[d][s];
        if (sumdano[d][s])
          anorcullis[s]    = anolofc[d][s]/sumdano[d][s];
      }
    }
    printonestat(d, "R Kraut", crkraut, arkraut, anorkraut);
    printonestat(d, "R Cullis", crcullis, arcullis, anorcullis);
    printonestat(d, "Phasing Power", cphasepower, aphasepower, anophasepower);
  }
}

void Bp3likelihood::printonestat(const unsigned d, const string stat, const vector<double> &cstat,
    const vector<double> &astat,
    const vector<double> &anostat) const
{

  printf(" $TABLE: Dataset %s  %s vs. Resolution-centric, acentric, overall, anom:\n",
      xtal.sf[d].name.c_str(), stat.c_str());
  if (xtal.sf[d].anomalous)
  {
    printf("$GRAPHS: %s vs. Res:N:4,7,9,11,13:\n", stat.c_str());
    printf("       : %s vs. STOL2:N:5,7,9,11,13:\n", stat.c_str());
  }
  else
  {
    printf("$GRAPHS: %s vs. STOL2:N:4,7,9,11:\n", stat.c_str());
    printf("       : %s vs. STOL2:N:5,7,8,11:\n", stat.c_str());
  }

  double totalcstat(ZERO);
  unsigned totalcnshl(0);

  double totalastat(ZERO);
  unsigned totalanshl(0);

  double totalanostat(ZERO);
  unsigned totalanonshl(0);

  printf("$$\n Bin   LoRes  HiRes    Res   STOL2    Refls   Centr   Refls  Acentr");
  printf("   Refls    All");

  if (xtal.sf[d].anomalous)
    printf("    Refls   Anom  $$\n$$\n");
  else
    printf("  $$\n$$\n");
  printf("\n");
  for (unsigned s   = 0; s < xtal.sf[d].nbins; s++)
  {
    double statbin  = ((anshl[d][s] + cnshl[d][s]) ?
        (astat[s]*anshl[d][s]+cstat[s]*cnshl[d][s])/
        ((double)anshl[d][s]+cnshl[d][s]) : ZERO);
    printf(" %2d  %6.2f  %6.2f %6.2f  %7.5f  %5d  %7.4f  %5d  %7.4f  %5d  %7.4f",
        s+1, xtal.sf[d].Getlores(s), xtal.sf[d].Gethires(s), xtal.sf[d].Getavres(s),
        xtal.sf[d].astolsq[s],
        cnshl[d][s],cstat[s],anshl[d][s],astat[s],anshl[d][s] + cnshl[d][s], statbin);
    if (xtal.sf[d].anomalous)
      printf("  %5d  %7.4f\n", sanonshl[d][s], anostat[s]);
    else
      printf("\n");

    totalcstat     += cstat[s]*cnshl[d][s];
    totalcnshl     += cnshl[d][s];
    totalastat     += astat[s]*anshl[d][s];
    totalanshl     += anshl[d][s];
    if (xtal.sf[d].anomalous)
    {
      totalanostat += anostat[s]*sanonshl[d][s];
      totalanonshl += sanonshl[d][s];
    }
  }
  double totcstat   = (totalcnshl) ? totalcstat/((double)totalcnshl) : ZERO;
  double totastat   = (totalanshl) ? totalastat/((double)totalanshl) : ZERO;
  double totstat    = ((totalanshl+totalcnshl) ?
      (totalastat+totalcstat)/((double)(totalanshl+totalcnshl))
      : ZERO);
  printf("$$\nTOTAL\n");  
  printf("     %6.2f  %6.2f                  %5d  %7.4f  %5d  %7.4f  %5d  %7.4f",
      xtal.sf[d].Getlores(0), xtal.sf[d].Gethires(xtal.sf[d].nbins-1), totalcnshl,
      totcstat, totalanshl, totastat, totalanshl+totalcnshl, totstat);

  if (xtal.sf[d].anomalous)
  {
    double totano   = (totalanonshl) ? totalanostat/totalanonshl : ZERO;
    printf("  %5d  %7.4f\n", totalanonshl, totano);
  }		 
  else
    printf("\n");

  printf("\n\n");
  fflush(stdout);  
}

void Bp3likelihood::print(const int cyc) const
{
  if (stats)
  {
    printfom();

    if ( (target     != "SAD") && (target != "MSRS") )
      printstat();
  }

  if ( (refinescale) && (target != "SAD") && (target != "MSRS") )
    xtal.printscale();

  if (refineluzzati)
    xtal.printLuzzati(target, cyc);

  mdl.print();
}
