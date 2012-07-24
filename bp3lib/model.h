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

#ifndef MODEL_H
#define MODEL_H	1

#include "misc.h"

class Scatter
{
  // Scatter holds normal scattering factors approximated by a Gaussian
  // distribution a[4], b[4], c - see $CLIBD/atomsf.lib for more info and
  // scattering factors dependent on wavelength (fp[w], fpp[w]).

  friend class Afrolikelihood;
  friend class Bp3likelihood;
  friend class Gcxlikelihood;
  friend class Multicomblikelihood;
  friend class Likelihood;
  
 public:
  Scatter(const string = "", const unsigned = 1, const unsigned = 1);
  ~Scatter(){;}

  Scatter &Setname(const string namein)
    {name = namein; return *this;}
  Scatter &Setnwave(const unsigned n)   { nwave = n; resize(); return *this; }
  Scatter &Setrefinefp (const unsigned, const bool);   
  Scatter &Setrefinefpp(const unsigned, const bool);
  Scatter &Setfp (const vector<double>);
  Scatter &Setfpp(const vector<double>);
  Scatter &Setfp (const unsigned, const double);
  Scatter &Setfpp(const unsigned, const double);
  
  double   Getfp (const unsigned w)      const {return fp[w];        } 
  double   Getfpp(const unsigned w)      const {return fpp[w];       }
  unsigned Getnwave()                    const {return nwave;        }  
  unsigned Getcrystal()                  const {return xtal;         }
  bool    Getrefinefp (const unsigned w) const {return refinefp[w];  } 
  bool    Getrefinefpp(const unsigned w) const {return refinefpp[w]; }  
  string  Getname() {return name;}  

  void print(const bool = false) const;

 private:
  string name;                        // atom name

  unsigned nwave;
  unsigned xtal;
  work a[4], b[4];                    // form factors a, b, and c
  work c;                             // which are just a function of name
  
  vector<double>  fp, fpp;            // fp, fpp and flags for refining
  vector<bool>  refinefp;             // them as a function of wavelength
  vector<bool>  refinefpp;

  Scatter &Setabc();                  // set forma, formb, formc and fp, fpp
                                      // given the atom name
  
  Scatter &resize();                  // resize vectors to appropriate
                                      // number of wavelengths
};

class Site
{
  friend class Afrolikelihood;
  friend class Bp3likelihood;
  friend class Gcxlikelihood;
  friend class Multicomblikelihood;
  friend class Likelihood;
  
 public:
  Site(const unsigned  = 0, const double = ZERO, const double = ZERO,
       const double = ZERO, const bool = true, const bool = true,
       const bool = true);
  ~Site(){;}

  Site &Setx(const unsigned i, const double value)
    { x[i] = between(MINX, value, MAXX); return *this; }
  Site &Setnumber (const unsigned  value)
  { number        = value; return *this;}
  Site &Setrefinex(const unsigned i, const bool value)
    { refinex[i]    = value; return *this;}
  double Getx(const unsigned i) const {return x[i]; }
  bool Getrefinex(const unsigned i) const {return refinex[i]; }
  vector<double> ortho(vector<vector<work> > &) const;
  void invert(const unsigned sgnumber);
  void print() const;

 private:
  unsigned number;
  double x[3];
  bool refinex[3];
  // fractional atomic coordinates
  double MINX, MAXX;
};

class Atom
{
  friend class Afrolikelihood;
  friend class Bp3likelihood;
  friend class Gcxlikelihood;
  friend class Multicomblikelihood;
  friend class Likelihood;
  
 public:
  Atom();
  
  Atom(const string, const unsigned, const unsigned,
       const unsigned, const bool, const vector<double> &,
       const double, const bool = true, const bool = true);
  ~Atom(){;}

  Atom &checkuaniso(const string);
  Atom &convertuaniso(const double);
  bool posdefu();

  Atom &Setname(const string namein)
  {name = namein; return *this;}

  Atom &Setuaniso(const double val, const unsigned i)
    {uaniso[i] = between(MINAB, val, MAXAB);return *this;}
  Atom &Setbiso(const double value)
    {biso = between(MINB, value, MAXB); return *this;}
  Atom &Setocc (const double value)
    {occ = between(MINO, value, MAXO); return *this;}
  Atom &Setnumber (const double value)
  {number = between(MINO, value, 500.0); refinenumb       = true  ; return *this;}  
  Atom &Setnumb   (const double value)
  {number = between(MINO, value, 500.0); refinenumb       = true  ; return *this;}  
  Atom &Setmino (const double value)
    {MINO = between(ZERO, value, MAXO); return *this;}
  Atom &Setminb (const double value)
    {MINB = between(ZERO, value, MAXB); return *this;}
  Atom &Setrefinebfac(const bool value)   { refinebfac    = value ; return *this;}
  Atom &Setrefinenumb (const bool value)  { refinenumb    = value ; return *this;}
  Atom &Setrefineocc (const bool value)   { refineocc     = value ; return *this;}

  string Getname()             const { return name;       }
  unsigned Getnsite()          const { return nsite;      }
  double Getbiso()             const { return biso;       }
  double Getuaniso(unsigned b) const { return uaniso[b];  }
  double Getocc()              const { return occ;        }
  double Getnumber()           const { return number;     }
  double Getnumb()             const { return number;     }
  unsigned Getcrystal()        const { return crystal;    }
  bool Getisotropic()          const { return isotropic;  }
  bool Getrefinebfac()         const { return refinebfac; }
  bool Getrefineocc()          const { return refineocc;  }
  bool Getrefinenumb()         const { return refinenumb; }

  void print() const;

 private:
  string name;
  unsigned nsite;           // index of the site of the atom
  double biso;
  double uaniso[6];         // uaniso[0] = U11, uaniso[1] = U12, uaniso[2] = U13
  double occ;               // uaniso[3] = U22, uaniso[4] = U23, uaniso[5] = U33
  unsigned crystal;         // the crystal the atom is in
  unsigned nform;           // index of the scattering factors for the atom 
  double number;

  bool isotropic;           // flag indicating whether we have isotropic
                            // or anisotropic bfactors
  bool refinebfac, refineocc, refinenumb;
  // maximum allowable values for parameters
  // occupancy parameters
  double MINO, MAXO;
  // bfactor parameters
  double MINB, MAXB;
  // anisotropic bfactor parameters
  double MINAB, MAXAB;
};

// Model class is just a class to hold all the model
// parameters (Site, Atom, and Scatter classes)

class Model
{
 public:
  Model(){;}
  ~Model(){;}

  vector<Atom> atom;
  vector<Scatter> form;
  vector<Site> site;

  void inverthand(const unsigned);
  
  void print(const bool = false) const;
};

#endif
