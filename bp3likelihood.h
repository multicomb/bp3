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

#ifndef BP3LIKELIHOOD_H
#define BP3LIKELIHOOD_H	1

#include "misc.h"
#include "cmtzlib.h"
#include "crystal.h"
#include "likelihood.h"
#include "matrix.h"
#include "model.h"

class Bp3likelihood : public Likelihood
{
 public:
    Bp3likelihood(Model &, Crystal &);
    ~Bp3likelihood(){;}
    void   setup(const bool = false, const bool = true, const bool = false);
    unsigned Setiterations() const;
    void   checkatoms();
    void   finitedifftest();                 // finite difference function test    
    double function(const bool = false);     // general function
    double grad(const bool = false, const bool = false); // general derivative
    
    virtual double gradient(const vector<double> &pars, vector<double> &grd)
    {double func(grad()); grd = Getgradient(); return func;}

    virtual double function1dim (const double stepsize, const vector<double> &pars,
				 const vector<double> &direction)
    {shift(stepsize,pars,direction); directsfcalc(false,updatesigmah); return function();}   

#if 0
    virtual double gradient1dim (double &dfunc, vector<double> &grd,
				 const double stepsize,
				 const vector<double> &pars,
				 const vector<double> &direction)
    { 
      std::string tag;

      tag = Tgrad1dim.start("shift");
      shift(stepsize,pars,direction); 
      Tgrad1dim.stop();

      tag = Tgrad1dim.start("directsfcalc");
      directsfcalc(false,updatesigmah); 
      Tgrad1dim.stop();

      tag = Tgrad1dim.start("grad");
      double func(grad());
      Tgrad1dim.stop();

      tag = Tgrad1dim.start("Getgradient");
      grd = Getgradient();
      Tgrad1dim.stop();

      tag = Tgrad1dim.start("dotprod");
      dfunc = dot_prod(grd,direction); 
      Tgrad1dim.stop();
      return func;
    }
#else
    double gradient1dim (double &dfunc, vector<double> &grd,
				 const double stepsize,
				 const vector<double> &pars,
				 const vector<double> &direction);
#endif

    virtual void print(const int = -1) const;
    void   otherhand();
    double norm() 
    {double nor(ZERO); vector<double> gra = Getgradient();
      for (unsigned i = 0; i < gra.size(); i++) nor += gra[i]*gra[i];
      return sqrt(nor);}
      
    // set functions for parser
    Bp3likelihood &Setshelthres(const unsigned input)
    {shelthres    = input; return *this;}
    Bp3likelihood &Setthreshold(const double input)
    {threshold    = input; return *this;}
    Bp3likelihood &Setsheldrick(const bool input)
    {sheldrick    = input; return *this;}
    Bp3likelihood &Setdifftol(const double input)
    {difftol      = input; return *this;}
    
    bool     Getsheldrick()    const {return sheldrick;    }
    unsigned Getmad()          const {return xtal.Getmad();}
  private:
    unsigned shelthres, cd;
    bool sheldrick;
    double threshold;                        // f/sigma threshold for 1D -> 2D switch
    double difftol;                          // finite difference test tolerance
    double minvar;  
  
    // stats
    vector<vector<double> >   aisolofc, cisolofc, anolofc;
    vector<vector<double> >   asumfph, csumfph, anosumfph;
    vector<vector<double> >   asumfh, csumfh, sumimfh;
    vector<vector<double> >   asumdiff, csumdiff, sumdano;
    vector<vector<unsigned> > sanonshl;

    double overallfom;
    void initstat();
    void printstat() const;
    void printonestat(const unsigned, const string, const vector<double> &,
		      const vector<double> &, const vector<double> &) const;
    
    // likelihood refinement functions
    double uncorrelatedfunction(const bool check); 
    double uncorrelatedgradient(const bool check, const bool mtz = false);
    double sadfunction(const bool check);
    double sadgradient(const bool check = false, const bool mtz = false);
    double sadgradient_gold(const bool, const bool);
    double madfunction(const bool check);
    double madgradient(const bool check = false, const bool mtz = false);

    // utility functions
    Bp3likelihood &Setgauss();
    void madphasingdefaults();    
};
#endif
