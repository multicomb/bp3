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

#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H	1

#include "fstream"
#include "misc.h"
#include "cmtzlib.h"
#include "crystal.h"
#include "matrix.h"
#include "model.h"
#include "tabfunc.h"

class Likelihood
{
  // methods needed for minimizer
 public:
  Likelihood(Model &, Crystal &);  // to use minimizer, the Likelihood class doesn't
  virtual ~Likelihood(){;}         // need Model or Crystal classes!

  // returns the maximum step to take in the line search - set to 1.0 if you
  // are unsure.

  double Getstepmax(const vector<double> &) const;

  // returns the (double) parameters that will be refined in a vector  

  vector<double> Setpars() const;  
  unsigned Getnpars() const {return npars; }

  // returns the function value and gradient vector of the parameters

  virtual double gradient (const vector<double> &, vector<double> &) = 0;
  
  // given a stepsize, the parameters and a search direction, it returns the function
  // values

  virtual double function1dim (const double,
			       const vector<double> &,
			       const vector<double> &)               = 0;
  
  // given a stepsize, the parameters and a search direction, it returns the function
  // value, the directional derivative (dfunc = grad dot_prod direction) and the
  // gradient vector (grad)

  virtual double gradient1dim (double &, vector<double> &,
			       const double,
			       const vector<double> &,
			       const vector<double> &d)              = 0;
  
  virtual void print(const int = -1) const                           = 0;
  void Setverbose(const unsigned input)  {verbose = input;   }
  unsigned Getverbose()  const           {return verbose;    }
  void Setwarn(const bool        input)  {warn    = input;   }
  bool Getwarn()  const                  {return warn;       }

  bool     Getnooutput() const           {return nooutput;   }
  double   Getfunction() const           {return likelihood; }
  double   Getbeta() const               {return beta; }
    
  vector<double> Getgradient() const;
  
 protected:
  double likelihood;
  // the likelhood variable holds the current function value in the likelihood class
  unsigned verbose;
  // verbose flag. 0 implies minimal messages, 1 implies more output, > 1 even more. 
  bool warn;  // if set to true, don't exit.

  // methods needed for crystallography
 public:
  // set functions for parser
  Likelihood &Setcommand(const char *input)
  {command       = input; return *this;}
  Likelihood &Setnosmall(const bool input)
  {nosmall       = input; return *this;}
  Likelihood &Setfill(const bool input)
  {fill          = input; return *this;}
  Likelihood &Setmaxocc(const double input)
  {maxoccshift   = input; return *this;}
  Likelihood &Setmaxluzz(const double input)
  {maxluzzshift  = input; return *this;}
  Likelihood &Setmaxbfac(const double input)
  {maxbfacshift  = input; return *this;}
  Likelihood &Setmaxxyz(const double input)
  {maxxyzshift   = input; return *this;}
  Likelihood &Setmaxx(const double input)
  {maxxshift     = input; return *this;}
  Likelihood &Setmaxy(const double input)
  {maxyshift     = input; return *this;}
  Likelihood &Setmaxz(const double input)
  {maxzshift     = input; return *this;}
  Likelihood &Setminocc(const double value)
  {for (unsigned i = 0; i < mdl.atom.size(); i++)
      mdl.atom[i].Setmino(value); return *this;}
  Likelihood &Setrefineluzzati(const bool input)
  {refineluzzati = input; return *this; }
  Likelihood &Setrefinescale(const bool input)
  {refinescale   = input; return *this; }
  Likelihood &Setrefineocc(const bool input)
  {refineocc     = input; return *this; }
  Likelihood &Setrefinenumb(const bool input)
  {refinenumb    = input; return *this; }
  Likelihood &Setrefinexyz(const bool input)
  {refinexyz     = input; return *this; }
  Likelihood &Setrefinebfac(const bool input)
  {refinebfac    = input; return *this; }
  Likelihood &Setrefineall(const bool input)
  {refineall     = input; return *this;}
  Likelihood &Setbeta(const double input)
  {beta          = input; return *this;}
  void applybeta() { //cout << "mlhl beta debug: " << dLddmod[0][0].size() << " "<< dLddmod[0][1].size() << endl;
    for (unsigned s = 0; s < xtal.sf[0].nbins; s++) {  if (xtal.pluz.size()) if (xtal.pluz[0].size()) xtal.pluz[0][s] *=beta; //cout << "here!! "  << endl;
      if (xtal.sf[0].dmod.size()) for (unsigned m                  = 0; m < xtal.sf[0].dmod[s].size(); m++) xtal.sf[0].dmod[s][m] *=beta;
    }
  }

  Likelihood &Setoutputhcalc(const bool input)
  {outputhcalc   = input; return *this;}
  Likelihood &Setmtzfb(const string input)
  {mtzfb         = input; return *this;}
  Likelihood &Setmtzpb(const string input)
  {mtzpb         = input; return *this;}
  Likelihood &Setmtzfdiff(const string input)
  {mtzfdiff      = input; return *this;}
  Likelihood &Setmtzpdiff(const string input)
  {mtzpdiff      = input; return *this;}
  Likelihood &Setmtzfcomb(const string input)
  {mtzfcomb      = input; return *this;}
  Likelihood &Setmtzpcomb(const string input)
  {mtzpcomb      = input; return *this;}
  Likelihood &Setmtzf(const string input)
  {mtzf          = input; return *this;}
  Likelihood &Setmtzsigf(const string input)
  {mtzsigf       = input; return *this;}
  Likelihood &Setmtzfa(const string input)
  {mtzfa         = input; return *this;}
  Likelihood &Setmtzsigfa(const string input)
  {mtzsigfa      = input; return *this;}
  Likelihood &Setmtzea(const string input)
  {mtzea         = input; return *this;}
  Likelihood &Setmtzsigea(const string input)
  {mtzsigea      = input; return *this;}
  Likelihood &Setmtzalpha(const string input)
  {mtzalpha      = input; return *this;}
  Likelihood &Setmtzfom(const string input)
  {mtzfom        = input; return *this;}
  Likelihood &Setmtzhla(const string input)
  {mtzhla        = input; return *this;}
  Likelihood &Setmtzhlb(const string input)
  {mtzhlb        = input; return *this;}
  Likelihood &Setmtzhlc(const string input)
  {mtzhlc        = input; return *this;}
  Likelihood &Setmtzhld(const string input)
  {mtzhld        = input; return *this;}
  
  void checkocc();
  void checkfrac(const vector<bool> &);
  void Setmaxshifts();
  void checkatoms();
  void Settitle(const string input)           {title         = input;   }
  void Settarget(const string input);
  void Setstats(const bool input)             {stats         = input;   }
  void Setallin(const bool input)             {allin         = input;   }
  void Setonly(const bool input)              {only          = input;   }
  void Setinterpolate(const bool input)       {interpolate   = input;   }
  void Setoutput(const string input)          {output        = input;   }
  void Setnooutput(const bool input)          {nooutput      = true;    }
  void Setmode(const string input)            {mode          = input;   }
  void Sethand(const unsigned input)          {hand          = input;   }
  void Setcweightsize(const unsigned input)   {cweight.resize(input);   }
  void Setpweightsize(const unsigned input)   {pweight.resize(input);   }
  void Setafweightsize(const unsigned input)  {afweight.resize(input);  }
  void Setsadweightsize(const unsigned input) {sadweight.resize(input); }
  void Setupdatesigmah(const bool input)      {updatesigmah  = input;   }
  
  string Gettitle()  const                    {return title;           }
  string Getmode()   const                    {return mode;            }
  string Gettarget() const                    {return target;          }
  unsigned Gethand() const                    {return hand;            }
  bool Getrefineall() const                   {return refineall;       }

  void writescript(const string = "bp3") const;
  void writepdb() const;
  void writexml() const;
  void outputintensities();
  void outputscafile() const;
  void outputhklfile() const;
  
 protected:
  Model &mdl;
  Crystal &xtal;
  TabFunc<double> tab;
  bool allin, only;                        // give input mtz columns in output
  bool interpolate, stats, nooutput, refineall, fill, outputhcalc;
  string command;
  string mode, protocol;
  string target;                           // target function requested
  string title;                            // title for MTZ file
  string output;
  unsigned hand;                           // run other hand (or not)
  unsigned npars;                          // number of parameters to refine
  // linking minimizer and likelihood function
  vector<string> type;
  void shift(const double, const vector<double> &, const vector<double> &);
  double xsshift(const double, const string &) const;
  double maxoccshift, maxluzzshift, maxxyzshift, maxbfacshift, maxuanoshift;
  double maxxshift, maxyshift, maxzshift;
  double maxfpshift, maxfppshift, maxkscaleshift, maxbscaleshift, maxnumbshift;
  double beta;
  
  void Settype();

  // gaussian integration variables
  vector<double> gsin, gcos, pweight;      // phase integral - acentric 1D
  vector<double> cnode, cweight;           // amplitude integral - centric
  vector<double> afnode, afweight;         // amplitude integral - acentric
  vector<double> sadsin, sadcos, sadweight;// sad integral

  // variables for matrix inversion
  vector<double> evectors, evalues, lawork;
  vector<int> iwork;
  int lwork, liwork, info;    
  bool hermitianinverse(Matrix &, Matrix &, Matrix &, Matrix &, double &, const bool = false);
  bool inverse(Matrix &, Matrix &, double &);
  void inverse2by2(Matrix &, Matrix &, double &);
  void hermitianmatrixprod(Matrix &, Matrix &,
			   const Matrix &, const Matrix &,
			   const double *, const double *);
  void matrixprod(Matrix &, const Matrix &, const double *);

  // structure factor calculation
  void directsfcalc(const bool, const bool);
  double fp(const unsigned, const unsigned, const unsigned);
  double calcanomcorrel() const;
  vector<vector<work> > hrot;
  vector<work> htran;
  bool onlyanocalc;
  // flags indicating parameter type(s) whose derivatives are required
  bool refineluzzati, refineocc, refinescale, refinexyz, refinebfac;
  bool refinefp, refinefpp, refinenumb, updatesigmah; 

  // derivatives
  vector<vector<vector<double> > > sumfpp, difffpp, difffp, sumfpfpp;

  // Note! Derivatives of structure factors are indexed [d][r] to save space if
  // there is no anomalous scattering.
  vector<vector<double> > dLdAp;           // Re(F+) - or Re(F) if F+=F=
  vector<vector<double> > dLdBp;           // Im(F+) - or Im(F) if F+=F-
  
  vector<vector<double> > dLdAm;           // Re(F-) (if anomalous)
  vector<vector<double> > dLdBm;           // Im(F-) (if anomalous)

  vector<vector<double> > dLdsigmah;       
  vector<vector<double> > dLdsumfpp;       
  vector<vector<double> > dLdsumfpfpp;       
  vector<double>          dLdocc;          // heavy atom parameters
  vector<double>          dLdnumb;
  vector<double>          dLdbiso;
  vector<vector<double> > dLduaniso;
  vector<vector<double> > dLdfp;
  vector<vector<double> > dLdfpp;
  vector<vector<double> > dLdx;

  vector<vector<double> > dLddluz;         // Luzzati parameters
  vector<vector<double> > dLdeluz;
  vector<vector<double> > dLdadluz;
  vector<vector<double> > dLdsdluz;
  vector<vector<double> > dLdpluz;
  vector<vector<vector<double> > > dLddmod;
   
  vector<double> dLdkscale;                // Scale parameters
  vector<double> dLdbscale;

  // SAD/MAD matrices
  vector<Matrix> covmodel, covinvmodel;
  Matrix redAdsd, redAdsigh, redAdsfpp, redAddmod;
  Matrix recov, recovinv;
    
  void resize(const bool outputmtz = false);
  
  // likelihood functions
  double multsirasgradient(const bool check = false, const bool mtz = false);
  double multsirasfunction(const bool check) { return multsirasgradient(); }
  double pavolsadgradient(const bool check = false, const bool mtz = false);
  double pavolsadfunction(const bool check) { return pavolsadgradient(); }
  double rejectionprobability;
  bool nosmall;
  // stats
  vector<double> afom, cfom;
  vector<vector<unsigned> > anshl, cnshl;
  void printfom() const;
  
  // mtz writing
  vector<CMtz::MTZCOL * > colin, colout;  
  void setupmtz(CMtz::MTZ *, CMtz::MTZ *);
  void storeinitialcolumns(CMtz::MTZ *, vector<float> &,
			   const unsigned , const unsigned );
  string mtzpb, mtzfb, mtzfom, mtzhla, mtzhlb, mtzhlc, mtzhld;
  string mtzf, mtzsigf;
  string mtzfdiff, mtzpdiff, mtzfcomb, mtzpcomb;
  string mtzfa, mtzsigfa, mtzea, mtzsigea, mtzalpha;
    
};
#endif  
