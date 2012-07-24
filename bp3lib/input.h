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

#ifndef INPUT_H
#define INPUT_H	1

#include "crystal.h"
#include "likelihood.h"
#include "model.h"
#include "minimizer.h"
#include "ccp4_parser.h"                  // Peter Briggs' parser in c libraries

class Input
{
 public:
  Input(Model &, Crystal &, Likelihood &, Minimizer &);
  ~Input(){;}
   
 protected:
  // references to  model, crystal, likelihood and minimizer classes
  Model &mdl;
  Crystal &xtal;
  Likelihood &like;
  Minimizer &mini;

  bool xtalkeywords();
  void xtalkeyword();      // XTAL KEYWORDS
  void cellkeyword();
  void centkeyword();
  void acenkeyword();
  void iampkeyword();
  void iintkeyword();
  void invekeyword();
  void scalkeyword();
  void pluzkeyword();

  bool datakeywords();
  void dnamkeyword();      // DATA KEYWORDS
  void mtzikeyword();
  void scaikeyword();
  void unkfkeyword();
  void colukeyword();
  void formkeyword();
  void resokeyword();
  void binskeyword();
  void kscalekeyword();
  void bscalekeyword();
  void isoekeyword();
  void anoekeyword();
  void corekeyword();
  void sdlukeyword();
  void exclkeyword();

  bool modelkeywords();
  void modlkeyword();      // MODEL KEYWORDS
  void sitekeyword();
  void atomkeyword();
  void xyzkeyword();
  void numbkeyword();
  void occukeyword();
  void bisokeyword();
  void uanokeyword();
  void norekeyword();

  bool minimizerkeywords();
  void cyclkeyword();      // MINIMIZER KEYWORDS
  void searkeyword();
  void linekeyword();
  void normkeyword();
  void tparkeyword();
  void wolfkeyword();

  bool likelihoodkeywords();
  void targkeyword();
  void fillkeyword(); 
  void nodekeyword(); 
  void outikeyword();
  void outhkeyword();
  void outskeyword();
  void intekeyword();
  void verbkeyword();
  void warnkeyword();
  void titlkeyword();
  void labokeyword();
  void statkeyword();
  void allikeyword();
  void outpkeyword();
  void nohakeyword();
  void noutkeyword();
  void noupkeyword();
  void refakeyword() const;
  void mocckeyword() const;
  void mluzkeyword() const;
  void mbfakeyword() const;
  void mxyzkeyword() const;
  void betakeyword();
 
  void check();
  void atomset();
  void readpdb(const char*);
  void resetatomin()
  { atomname = ""; setsite = setb = seto = false; number = ONE; refo = refb = refnumb = true;}
  
  bool atomready()
  {  return (setsite) && (seto) && (atomname != "");}

  CCP4::CCP4PARSERARRAY *parser;    // Peter's parser/tokenizer structure

  vector<bool> frac;                // orthogonal or fractional coordinates given

  int nsite;                        // current site being parsed
  int natom;                        // current atom being parsed
  int ndata;                        // current data set being parsed
  int nxtal;                        // current crystal number being parsed
  int nwave;                        // current wavelength being parsed
  
  string atomname;
  vector<string> unknownfiletype;
  double occ;
  vector<double> bfac;
  bool inputamplitudes, inputintensities;
  bool norefocc, norefxyz, norefbfac;
  bool crystal, dataset;
  bool biso, bfacset;
  bool setsite, setwithsite, setb, seto, setnumb, setwithmodl, setwithatom;
  bool refo, refb, refnumb;
  double number;
  unsigned mbin, refbin;
};

#endif
