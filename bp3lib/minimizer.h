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

#ifndef MINIMIZER_H
#define MINIMIZER_H	1

#include "misc.h"
#include "likelihood.h"

class Minimizer
{
 public:
  Minimizer(Likelihood &like) : likelihood(like)
  { maxit  = 500; minit = 25; maxlineit = 15; search = "BFGS"; iter = reinit = 0; normtol = 10.0;
    partol = 0.0001; alpha  = 0.0001; beta = 0.975; line = "MORE"; stepmin = 1.0e-20; step = 0.00001;
    setit = false;}
  ~Minimizer(){;}
  void refine(const bool savehessian = false);

  unsigned  &Getiterations(){return maxit;}
  bool      &Getsetit()     {return setit;}

  Minimizer &Setiterations(const unsigned input){maxit   = input; setit = true; return *this;}
  Minimizer &Setnormtol   (const double input)
  {normtol = std::max(input, EPSILON); return *this;}
  Minimizer &Setstep     (const double input)
  {step = std::max(input, 0.00001); return *this;}
  Minimizer &Setpartol    (const double input)
  {partol  = std::max(input, EPSILON); return *this;}
  Minimizer &Setsearch    (const string input)
  { if (input == "BFGS" || input == "DFP") line    = input;
    else                                   line    = "BFGS";
    return *this;}
  Minimizer &Setline      (const string input)
  { if (input == "MORE" || input == "BRENT") line    = input;
    else                                     line    = "MORE";
    return *this;}
  Minimizer &Setalpha     (const double input)
  {alpha   = std::max(DSMALL,std::min(0.1, input)); return *this;}
  Minimizer &Setbeta      (const double input)
  {beta    = std::min(ONE-DSMALL,std::max(0.9,input)); return *this;}

 private:
  Likelihood &likelihood;
  bool setit;
  unsigned iter;                    // current iteration
  unsigned reinit;
  unsigned maxit;                   // maximum number of iterations
  unsigned minit;                   // minimum number of iterations
  unsigned maxlineit;               // maximum line search number of iterations
  unsigned nfe;                     // number of function evaluations
  string search;                    // choices are "BFGS" or "DFP"
  double firstfunc;                 // first function value
  double func, dfunc;               // function and directional derivative of nth it.
  double newfunc, dnewfunc;         // function and directional derivative of n+1th it.
  vector<double> grad, oldgrad;     // gradient  vector of nth and n-1th it.
  vector<double> pars, oldpars;     // parameter vector of nth and n-1th it.
  vector<double> direction;         // direction of search calculated
  vector<vector<double> > hessian;
  unsigned negativegamma;
  double norm, normtol, partol;
  void reinitializehessian()
  { Bp3Warning("Minimizer::reinitializehessian","reinitializing hessian");
    for(unsigned i     = 0; i < grad.size(); i++)
      {for (unsigned j  = 0; j < grad.size(); j++)
	  hessian[i][j] = ZERO;
	hessian[i][i]   = ONE;
      }}
  // line minimizer variables
  bool samegrad;                    // gradient calculated in lineminizer?
  bool steepest;
  bool bracket;                     // is a line minimum bracketted?
  unsigned lineconverge;
  string line;                      // choices are MORE and BRENT
  double alpha, beta;               // armijo/goldstein constants 
  double step, stepmin, stepmax;          // minimum and maximum for the line search
  double fbest, dfbest, stepbest;   // best parameters from line search
  // Search direction function(s)
  void bfgs();
  void dfp();
  // Line minimizer functions
  void lineminimizer();
  bool more();
  void choosestep(double &, double &, double &, double &, double &, 
		  double &, double &, double &, const double , const double );
  double brent(double &);
  void checkshift();
  // utility functions
  void resize();
  void print() const;  
};

#endif
