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

#ifndef CRYSTAL_H
#define CRYSTAL_H	1

#include "misc.h"
#include "sfdata.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>

class Cell
{
  friend class Afrolikelihood;
  friend class Bp3likelihood;
  friend class Bricslikelihood;
  friend class Gcxlikelihood;
  friend class Multicomblikelihood;
  friend class Crystal;
  friend class Likelihood;
 public:
  Cell(){name = ""; a = b = c = 0.0; alpha = beta = gamma = 90.0;}
  ~Cell(){;}
  
 private:
  // name of crystal
  string name;
  
  // cell edges and angles
  work a, b, c, alpha, beta, gamma;

  // reciprocal cell edges
  work astar, bstar, cstar;

  // volume
  work volume;
  
  // Orthogonal to fractional matrix
  vector<vector<work> > or2frac;
  // Fractional to orthogonal matrix
  vector<vector<work> > frac2or;

  // sin(theta) over lambda squared
  vector<work> stolsq;

  // Set functions for cell parameters (making sure they are positive)
 public:
  Cell &Setname(const string value){name = value; return *this;}
  Cell &Setcell(const float *);
  Cell &Seta(const work value, const work low = WEPSILON){
    a = (value > low)? value : low; return *this; }
  Cell &Setb(const work value, const work low = WEPSILON){
    b = (value > low)? value : low; return *this; }
  Cell &Setc(const work value, const work low = WEPSILON){
    c = (value > low)? value : low; return *this; }

  Cell &Setalpha(const work value, const work low = WEPSILON){
    alpha = (value > low)? value : low; return *this; }
  Cell  &Setbeta(const work value, const work low = WEPSILON){
    beta  = (value > low)? value : low; return *this; }
  Cell &Setgamma(const work value, const work low = WEPSILON){
    gamma = (value > low)? value : low; return *this; }

  string Getname() const {return name;}

 private:
  // void Checkcell(const Spacegroup &);
  void Setor2frac();
  void Setabcstar();
  void Setstolsq(const vector<vector<int> > &);
  void Setstolsq(const vector<vector<int> > &, const unsigned);
  void print(const bool = false) const;
};

class Spacegroup
{
  friend class Afrolikelihood;
  friend class Bp3likelihood;
  friend class Bricslikelihood;
  friend class Gcxlikelihood;
  friend class Multicomblikelihood;
  friend class Crystal;
  friend class Likelihood;

 public:
  Spacegroup(const int, const char* = NULL);
  Spacegroup(){;}
  ~Spacegroup(){;}
  string Getsystem()              {return system;}
  string Gethmname()              {return hmname;}
  string Getpointgroup()          {return pointgroup;}
  bool Getpolar(unsigned i) const {return polar[i];}
  bool specialpos(double &, double &, double &, double &) const;
  unsigned enantiomorph()   const;
  unsigned Getnumber()      const {return number;} 

 private:
  // Number of all and primitive symmetry elements
  unsigned NSYM, NSYMP;
  // spacegroup (Herman Mauguin) name
  string hmname;
  // Laue point group
  string pointgroup;
  // Crystal system
  string system;
  // CCP4 spacegroup number
  unsigned number;
  // Lattice type (P, A, B, C, F, I, R, H)
  char lattice;
  // polar axis?
  bool polar[3];
  // Symmetry matrix (transposed)
  vector<vector<vector<int> > > symrot;
  vector<vector<double> > symtran;

  void print(const bool = false) const;
  unsigned naniso() const;
  void checklattice() const;
  Spacegroup &Setsystem();
  Spacegroup &Setpolar();
};

class Crystal
{
  friend class Afrolikelihood;
  friend class Bp3likelihood;
  friend class Bricslikelihood;
  friend class Gcxlikelihood;
  friend class Multicomblikelihood;
  friend class Input;
  friend class Likelihood;
 
  public:
  // functions
    Crystal();
    ~Crystal(){;}

    Crystal &setwithmtz(const string type = "REFINE");
    Crystal &setwithunknownfiletype(const string type, const vector<string> name);
    Crystal &setwithmanymtz(const string type = "REFINE", const bool amp = true);
    Crystal &setwithmanysca(const string type = "REFINE");
    Spacegroup sg;                                           // spacegroup info
    vector<Sfdata> sf;                                       // structure factor data!

    // clipper objects
    Crystal &Setclipper (const bool input ){clipper  = input; return *this;}
    Crystal &setwithpdb(char * filename);
    Crystal &setwithpdb_ref(char * filename);
    //Crystal &Setmrqdata(int ndata);
    void Setspg_trial(const int input)          {spg_trial    = input;   }
    clipper::MiniMol model_base, model_result, ref_model;
    bool phasematch;
	int spg_trial;
    
    unsigned bin(const unsigned d, const unsigned r) const   // bin of dataset/ref
    {return between((unsigned)0, (unsigned) (work)((cell[sf[d].nxtal].stolsq[r]  -
						    sf[d].lowsq)*sf[d].vshell), sf[d].nbins-1);}
    unsigned Getndata() const {return sf.size();}
    bool     Getmad() const {return mad;}
    bool     Getinvert() const {return invert;}
    string Getname(unsigned c) const {return cell[c].Getname();}
    vector<float> Getcell(const unsigned d) const;
    void print() const;
    void printLuzzati(const string targ = "", const int = -1) const;
    void printscale() const; 
    double Calculatebfactor(const unsigned, const unsigned, const unsigned, const unsigned = 0) ;
    void Calculatematthews(const unsigned, const unsigned, double, unsigned*, double *) const;

    // stuff pavol's functions
    Crystal &Setpluz(const unsigned p, const vector<double> &input) {pluz.resize(p+1); for (unsigned i = 0; i < input.size(); i++)
        pluz[p].push_back(between(0.0,input[i],25.0)); return *this;}
    Crystal &Setpluz(const unsigned p, const unsigned s, const double input)
                                                                    {pluz[p][s]     = between(0.0,input,10.0); return *this;}
    Crystal &Setdefaultpluz(const unsigned p, const double input)   {defaultpluz[p] = between(0.0,input,10.0); return *this;}
    Crystal &Setrefinep(const unsigned p, const bool input)         {refinep.resize(p+1); refinep[p]  = input; return *this;}
    Crystal &Setuserpluz(const unsigned p, const bool input)        {userpluz[p]    = input; return *this;}
    double  Getpluz(const unsigned p, const unsigned s)             {return pluz[p][s];}
    vector<double>  Getpluz(const unsigned p)                       {return pluz[p];}

    // Target's needed by parser
    Crystal &Setverbose(const unsigned input)   {verbose       = input; return *this;}
    Crystal &Setwarn(const bool input)          {warn          = input; return *this;}
    Crystal &Setrscale  (const unsigned input)  {rscale        = input; return *this;}
    Crystal &Setonlycentrics (const bool input ){onlycentrics  = input; return *this;}
    Crystal &Setonlyacentrics(const bool input ){onlyacentrics = input; return *this;}
    Crystal &Setinvert(const bool input )       {invert        = input; return *this;}

 private:
    vector<Cell> cell;                                 // Cell info
    vector<string> pdbfilename;
    unsigned maxref;                                   // max number of reflections
    unsigned maxselref;                                // max index of selected refls
    unsigned verbose;                                  // verbose flag
    bool warn;                                         // warn flag
    vector<vector<int> >  miller;                      // Miller Indices
    vector<bool> selected;                             // reflections in refinement
    vector<work> centricphase;                         // centric phase restriction
    vector<unsigned> epsilon;                          // epsilon
    vector<bool> centric;                              // centric reflection ?
    bool rscale;
    bool onlyacentrics, onlycentrics;                  // flags to only refine
                                                       // acentrics/centrics
    bool clipper;                                      // using clipper's data structures?
    bool mad, heavyref, invert;
    int nd;
    // give the dataset corresponding to a given crystal and wavelength
    vector<vector<unsigned> > dataset;
    
    // pavol's functions
    vector<vector<double> > pluz;
    vector<bool> refinep, userpluz;
    vector<double> defaultpluz;
    
    // utility functions
    void resize();
    work Getres(const unsigned d, const unsigned r) const
    {return ONE/(TWO*sqrt(cell[sf[d].nxtal].stolsq[r]));}
    void binweights(const unsigned, const unsigned, unsigned &,
		    unsigned &, double &, double &) const;
    bool    Sysabs(const unsigned) const;
    Crystal &Setepsilon(const unsigned);               // set functions
    Crystal &Setcentric(const unsigned);
    Crystal &Setreso();
    Crystal &Setdataset();
    Crystal &Setselected(const string = "REFINE");
    Crystal &Setreflectionsphased();
    Crystal &Setnshl(const unsigned r);
    Crystal &Setnormalization();
    Crystal &SetLuzzati();
    Crystal &Guessdluz();                              // initial guess for luzzati
    Crystal &Guesseluz();                              // parameters
    Crystal &Guessdmod();
    Crystal &Guessadluz();
    Crystal &Guesssdluz();
    Crystal &Guesssdluzref();
};

#endif
