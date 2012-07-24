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

#ifndef SFDATA_H
#define SFDATA_H	1

#include <math.h>
#include "misc.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>

class Sfdata
{
  friend class Afrolikelihood;
  friend class Bp3likelihood;
  friend class Bricslikelihood;
  friend class Gcxlikelihood;
  friend class Multicomblikelihood;
  friend class Crystal;
  friend class Likelihood;
 public:
  Sfdata(const unsigned = 25, const unsigned = 2000);
  ~Sfdata(){;}
  Sfdata &Setscain (const string   val) {
    if (mtzin.size()) Bp3Error("Sfdata::Setscain","Can not set both input sca and mtz");
    else scain  = val; return *this;}
  Sfdata &Setmtzin (const string   val) {
    if (scain.size()) Bp3Error("Sfdata::Setmtzin","Can not set both input mtz and sca");
    else mtzin  = val; return *this;}
  Sfdata &Setname  (const string   val) {name   = val; return *this;}
  Sfdata &Settype  (const string   val) {type   = val; return *this;}
  Sfdata &Setnxtal (const unsigned val) {nxtal  = val; return *this;}
  Sfdata &Setnwave (const unsigned val) {nwave  = val; return *this;}

  Sfdata &Setsdata(const string val) {sdata     = val; return *this;}
  Sfdata &Setsdev (const string val) {sdev      = val; return *this;}
  Sfdata &Setsdatap(const string val) {sdatap   = val; return *this;}
  Sfdata &Setsdevp (const string val) {sdevp    = val; return *this;}
  Sfdata &Setsdatam(const string val) {sdatam   = val; return *this;}
  Sfdata &Setsdevm (const string val) {sdevm    = val; return *this;}
  Sfdata &Setsfmodel(const unsigned m, const string val)
  {sfmodel.resize(m); sfmodel[m-1]              = val; return *this;}
  Sfdata &Setspmodel(const unsigned m, const string val)
  {spmodel.resize(m); spmodel[m-1]              = val; return *this;}
  Sfdata &Setsfom(const string val) {sfom       = val; return *this;}
  Sfdata &Setsphi(const string val) {sphi       = val; return *this;}
  Sfdata &Setshla(const string val) {shla       = val; return *this;}
  Sfdata &Setshlb(const string val) {shlb       = val; return *this;}
  Sfdata &Setshlc(const string val) {shlc       = val; return *this;}
  Sfdata &Setshld(const string val) {shld       = val; return *this;}

  Sfdata &Sethires    (const work val) {hires   = val; return *this;}
  Sfdata &Setlowres   (const work val) {lowres  = val; return *this;}

  Sfdata &Setusernbins(const unsigned val)
  {if (val      > MAXBINS)
      MAXBINS   = val;
    nbins       = between(MINBINS, val, MAXBINS); return *this;}

  Sfdata &Setdefaultsdluz(const double val)
  {defsdluz     = between(MINSDLUZ, val, MAXSDLUZ); return *this;}

  Sfdata &Setkscale (const double val)
  {kscale       = between(MINK, val, MAXK); return *this;}
  Sfdata &Setbiso   (const double val)
  {biso         = between(MINBSCLISO, val, MAXBSCLISO); return *this;}

  Sfdata &Setdluz (const unsigned, const double);
  Sfdata &Seteluz (const unsigned, const double);
  Sfdata &Setadluz(const unsigned, const double);
  Sfdata &Setsdluz(const unsigned, const double);
  Sfdata &Setdmod (const unsigned, const unsigned, const double);

  Sfdata &Setdluz (const vector<double> &);
  Sfdata &Seteluz (const vector<double> &);
  Sfdata &Setadluz(const vector<double> &);
  Sfdata &Setsdluz(const vector<double> &);

  string Getname()        const {return name;}  
  string Gettype()        const {return type;}  
  string Getsdata()       const {return sdata;}
  string Getsdatap()      const {return sdatap;}
  string Getsdatam()      const {return sdatam;}
  string Getsdev()        const {return sdev;}
  string Getsdevp()       const {return sdevp;}
  string Getsdevm()       const {return sdevm;}
  string Getsfom()        const {return sfom;}
  string Getsphi()        const {return sphi;}
  string Getshla()        const {return shla;}
  string Getshlb()        const {return shlb;}
  string Getshlc()        const {return shlc;}
  string Getshld()        const {return shld;}
  string Getmtzin()       const {return mtzin;}
  string Getscain()       const {return scain;}

  unsigned Getnbins()     const {return nbins;}  
  unsigned GetMINBINS()   const {return MINBINS;}  
  bool     Getanomalous() const {return anomalous;}  

  double Getkscale()      const {return kscale;}  
  double Getbiso()        const {return biso;}      

  work datamean(const unsigned r) const
  {return ( (anomalous) ?
	    (anouse(r) ? HALF*(datap[r] + datam[r]) :
	     ( usem(r) ? datam[r] : datap[r]) ) :
	    datap[r] );}

  work hmean(const unsigned r) const
  {return ( (anomalous) ?
	    HALF*(fcalcp[r] + fcalcm[r]) :
	    fcalcp[r] );}

  work devmean(const unsigned r) const
  {return ( (anomalous) ?
	    (anouse(r) ? HALF*(devp[r] + devm[r]) :
	     ( usem(r) ? devm[r] : devp[r]) ) :
	    devp[r] );}

  work dano(const unsigned r)    const {return datap[r] - datam[r];}
  work sigdano(const unsigned r) const
  {return sqrt(devp[r]*devp[r]+devm[r]*devm[r]);}
  bool use (const unsigned r) const
  {return (anomalous) ? ( (datap[r] > NOTUSED) || (datam[r] > NOTUSED) ) :
      (datap[r]     > NOTUSED);}
  bool usem (const unsigned r) const
  {return (datam[r] > NOTUSED);}
  bool usep (const unsigned r) const
  {return (datap[r] > NOTUSED);}
  
  bool anouse (const unsigned r) const
  {return (datap[r] > NOTUSED) && (datam[r] > NOTUSED);}
      
  double Getdluz(const unsigned s)                    const {return dluz[s];}  
  double Geteluz(const unsigned s)                    const {return eluz[s];}      
  double Getadluz(const unsigned s)                   const {return adluz[s];}  
  double Getsdluz(const unsigned s)                   const {return sdluz[s];}  
  double Getdmod(const unsigned s, const unsigned m)  const {return dmod[s][m];}  

  Sfdata &Setrefinek (const bool val)   {refinek        = val; return *this;}
  Sfdata &Setrefineb (const bool val)   {refineb        = val; return *this;}
 
  Sfdata &Setrefined (const bool val)   {refined        = val; return *this;}
  Sfdata &Setrefinee (const bool val)   {refinee        = val; return *this;}
  Sfdata &Setrefinead(const bool val)   {refinead       = val; return *this;}
  Sfdata &Setrefinesd(const bool val)   {refinesd       = val; return *this;}
  Sfdata &Setrefinedmod(const bool val) {refinedmod     = val; return *this;}

  Sfdata &Setsigmacut(const double val)    {sigmacut    = val; return *this;}
  Sfdata &Setanosigmacut(const double val) {anosigmacut = val; return *this;}
  Sfdata &Setisosigmacut(const double val) {isosigmacut = val; return *this;}
  Sfdata &Setanormscut(const double val)   {anormscut   = val; return *this;}

  Sfdata &Setnbins(const unsigned input)
    {nbins = between(MINBINS, input , MAXBINS); return *this;}

  work Gethires(const unsigned s) const
  {return ONE/(TWO*sqrt((work)(s+1)/vshell + lowsq));} 
  work Getlores(const unsigned s) const
  {return ONE/(TWO*sqrt((work)(s  )/vshell + lowsq));}
  work Getavres(const unsigned s) const
  {return ONE/(TWO*sqrt(astolsq[s]));}
  
 private:
  // individual mtz or sca file for the data set (if given)
  string mtzin, scain;
  // user given name of the data set
  string name;
  // amplitude or intensity
  string type;
  // name of diffraction data in mtz file
  string sdata, sdev, sdatap, sdevp, sdatam, sdevm;
  vector<string> sfmodel, spmodel;
  string sfom, sphi, shla, shlb, shlc, shld;
  bool anomalous, noniso;
  // the diffraction data (a function of MAXREF)
  vector<work> datap, devp;
  vector<work> datam, devm;
  // Fcalc+/-, Pcalc+/-
  vector<double> fcalcp, pcalcp;
  vector<double> fcalcm, pcalcm;
  vector<vector<double> > fmodel, pmodel;
  vector<double> fom, phib, hla, hlb, hlc, hld;
  // clipper data sets
  clipper::HKL_info hkl_list, hkl_list_data, hkl_list_data_mtz, hkl_list_P1;
  // define objects
  clipper::HKL_data<clipper::data32::F_phi> Fna, Fref;
  clipper::HKL_data<clipper::data32::F_sigF> Fobs, Fobs_mtz, Fobs_P1;
  // sigma cut-off
  double sigmacut, anosigmacut, isosigmacut, anormscut;
  // number of reflections in this dataset
  unsigned nref;                        
  unsigned nxtal;                       // crystal dataset is in
  unsigned nwave;                       // wavelength # of dataset
  // resolution limits                  
  work hires, lowres;
  unsigned nbins;                       // number of bins/data set
  work lowsq, hisq, vshell;             // binning variables
  
  // normalization variables
  vector<double> sigman, sigmanref;
  vector<double> sigmah, sigmap, sigmadano;
  vector<unsigned> nshl, anonshl;
  vector<double> astolsq;
  // parameters to be refined
  double kscale, biso, baniso[6];
  vector<vector<double> > dmod;         // luzzati d parameter for models
  vector<double> dluz, eluz;            // luzzati (d and e)
  vector<double> adluz;                 // luzzati (anomalous d)
  vector<double> sdluz;                 // luzzati for sad data sets
  bool refinek, refineb, refinedmod;
  bool refined, refinee, refinead, refinesd, refinesdref;
  bool aniso;

  // statistics
  vector<double> fovers, danoovers;

  double defsdluz;
  
  // maximum allowable values for parameters
  // scale parameters
  double MINK, MAXK, MINBSCLISO, MAXBSCLISO;

  // luzatti parameters
  double MINDLUZ, MAXDLUZ, MINADLUZ, MAXADLUZ;
  double MINELUZ, MAXELUZ, MINSDLUZ, MAXSDLUZ;

  // bins
  unsigned MINBINS, MAXBINS, REFBINS;
  // resolution
  work MINRES, MAXRES;
  
  // utility functions
  void resize(const unsigned, const bool);
  void resizeclipper();
  void resizebin();
  double foversigma(const unsigned r) const {return (double)(datamean(r)/devmean(r));}
  double var(const double k, const unsigned s, const bool heavyref = false) const
  {return ( (heavyref) ?
	    (k*k*sigman[s] - dluz[s]*dluz[s]*sigmanref[s]) : 
	    (k*k*sigman[s] - dluz[s]*dluz[s]*sigmanref[s] - sigmah[s]) );}
  
  Sfdata &Setnbins()
    {nbins      = between(MINBINS, nref/REFBINS, MAXBINS); return *this;}
  Sfdata &Setshell();
  Sfdata &Setnref();
  void checkparameters();
  void printLuzzati(const string, const int = -1, const bool = false) const;
  void printsignal() const;
  void printnorm(const bool verbose = false) const;
  void print() const;
};

#endif
