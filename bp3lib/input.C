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

#include <ctype.h>
#include <fstream>
#include <iostream>
#include "input.h"

Input::Input(Model &model, Crystal &xtl, Likelihood &likelihood,
		   Minimizer &minimizer)
  : mdl(model), xtal(xtl), like(likelihood), mini(minimizer)
{
  // initialize variables before parsing
  crystal                 = dataset       = false;
  natom                   = ndata = nwave = nsite = nxtal = -1;
  // bin calculation
  mbin                    = 25;
  refbin                  = 2000;
  // atomic parameters
  setwithsite             = setwithmodl   = setwithatom   = false ;
  setsite                 = setb = seto   = setnumb       = false;
  norefocc                = norefxyz      = norefbfac     = false;
  atomname                = "";
  occ                     = 0.25;
  biso                    = true;
  bfac.resize(1, 25.0);

  inputamplitudes         = inputintensities              = false;
  
  // initialize ccp4 parser structure to store 50 tokens
  parser                  = (CCP4::CCP4PARSERARRAY *) CCP4::ccp4_parse_start(50);
}

void Input::check() 
{
  if (atomready())
    atomset();

   // checks if input is okay
  
  if (xtal.sf.size()         <= 0)
    Bp3Error("Input::check", "no data sets have been defined");

  /*
  if ( (xtal.sf.size()       == 1) && (xtal.sf[0].Getsdatam() == "") )
    Bp3Error("Input::check", "only one nonanomalous data set is given");
  */

  string inputtype            = xtal.sf[0].Gettype();
  for (unsigned d             = 1; d < xtal.sf.size(); d++)
    if (xtal.sf[d].Gettype() != inputtype)
      Bp3Error("Input::check", "input either all intensities or amplitudes - not both");
  
  bool onemtz(true);
  for (unsigned d             = 0; d < xtal.sf.size(); d++)
  {
    if (xtal.sf[d].Getmtzin().size() || xtal.sf[d].Getscain().size() || unknownfiletype.size())
      onemtz                  = false;

    if (!xtal.sf[d].Getsdatam().size())
    {
      if (xtal.sf[d].Getsdata().size() && !xtal.sf[d].Getsdatap().size())
	xtal.sf[d].Setsdatap(xtal.sf[d].Getsdata());
      if (xtal.sf[d].Getsdev().size() && !xtal.sf[d].Getsdevp().size())
	xtal.sf[d].Setsdevp(xtal.sf[d].Getsdev());
    }
  }
  
  
  if ((getenv("HKLIN")       == NULL) && onemtz)
    Bp3Error("Input::check", "mtz file not defined");
  if ((getenv("HKLIN")       != NULL) && !onemtz)
    Bp3Error("Input::check", "define one overall mtz file OR individual dataset(s) mtz file(s) - not both");

  if (inputamplitudes && inputintensities)
    Bp3Error("Input::check", "do not use both IAMP and IINT keywords!");

  if (inputamplitudes)
    for (unsigned d    = 0; d < xtal.sf.size(); d++)
      xtal.sf[d].Settype("AMPLITUDE");
  else if (inputintensities)
    for (unsigned d    = 0; d < xtal.sf.size(); d++)
      xtal.sf[d].Settype("INTENSITY");
    
  if (getenv("HKLIN") != NULL)
    xtal.setwithmtz("ALL");
  else
    for (unsigned d    = 0; d < xtal.sf.size(); d++)
      if (xtal.sf[d].Getmtzin().size())
      {
	  
        xtal.setwithmanymtz("ALL",false);
        break;
      }
      else if (xtal.sf[d].Getscain().size())
      {
        xtal.setwithmanysca("ALL");
        break;
      }
      else if (unknownfiletype[d].size())
      {
        xtal.setwithunknownfiletype("ALL",unknownfiletype);
        break;
      }  

  // some checking of input ie. should coordinates not be refined due
  // to special positions or polar spacegroups.

  /*
  if ( (like.Getmode()   == "PHASE") || (like.Getmode()   == "CHECK") )
    for (unsigned d       = 0; d < xtal.sf.size(); d++)
      xtal.sf[d].Setdefaultsdluz(0.65);
   */

  if (mdl.site.size() && mdl.atom.size())
  {
    // like.checkocc();

    printf("Making sure atoms have fractional coordinates less than ONE\n\n");
    like.checkfrac(frac);
//    like.Setmaxshifts();

    // fix x, y, and/or z for one atom if spacegroup is polar
    for (unsigned i         = 0; i < 3; i++)
      if (xtal.sg.Getpolar(i))
	if (mdl.site.size() > 0)
	  mdl.site[0].Setrefinex(i,false);

    for (unsigned i         = 0; i < mdl.atom.size(); i++)
    {
      unsigned s(mdl.atom[i].Getnsite());

      double x[3]           = {mdl.site[s].Getx(0),
			       mdl.site[s].Getx(1),
			       mdl.site[s].Getx(2)};

      if (norefxyz)
	mdl.site[s].Setrefinex(0,false).Setrefinex(1,false).Setrefinex(2,false);

      if (norefocc)
	mdl.atom[i].Setrefineocc(false);

      if (norefbfac)
	mdl.atom[i].Setrefinebfac(false);

      double occ            = mdl.atom[i].Getocc();

      if (xtal.sg.specialpos(x[0], x[1], x[2], occ))
      { 
	for (unsigned c     = 0; c < 3; c++)
	  mdl.site[s].Setx(c, x[c]);
	mdl.atom[i].Setocc(occ);
	printf("Warning: Site %d is on a special position\n", i+1);
	printf("Its occupancy is being divided by the multiplicity.\n\n");

      }
    }
  }
}

// XTAL keywords

bool Input::xtalkeywords()
{
  if (CCP4::ccp4_keymatch("XTAL", parser->keyword))
    xtalkeyword();
  else if (CCP4::ccp4_keymatch("CELL", parser->keyword))
    cellkeyword();
  else if (CCP4::ccp4_keymatch("CENT", parser->keyword))            
    centkeyword();
  else if (CCP4::ccp4_keymatch("ACEN", parser->keyword))            
     acenkeyword();
  else if (CCP4::ccp4_keymatch("IAMP", parser->keyword))            
     iampkeyword();
  else if (CCP4::ccp4_keymatch("IINT", parser->keyword))            
     iintkeyword();
  else if (CCP4::ccp4_keymatch("INVE", parser->keyword))            
     invekeyword();
  else if (CCP4::ccp4_keymatch("SCAL", parser->keyword))            
     scalkeyword();
  else if (CCP4::ccp4_keymatch("PLUZ", parser->keyword))            
    pluzkeyword();
  else
    return false;

  return true;
}

void Input::xtalkeyword()
{
  if (atomready())
    atomset();
  
  // update the crystal and wave counters
  nxtal++;
  crystal             = true;
  dataset             = false;
  nwave               = -1;
  resetatomin();

  // increment the cell for this crystal
  xtal.cell.push_back(Cell());
  xtal.pdbfilename.push_back("");
  frac.push_back(true);

  // check to see if we have no more than two tokens
  if (parser->ntokens > 2)
    Bp3Warning("Input::xtalkeyword", "Extra text in XTAL keyword will be ignored");

  if (parser->token[1].isstring)
    xtal.cell[nxtal].Setname(parser->token[1].fullstring);
  else
    Bp3Error("Input::xtalkeyword", "A string is expected for the XTAL's name");
}

void Input::cellkeyword()
{
  if (!crystal)
    Bp3Error("Input::cellkeyword","XTAL keyword must be given before CELL");
  
  // check to see if we have no more than 7 tokens
  if (parser->ntokens  > 7)
    Bp3Error("Input::cellkeyword", "Extra text in Cell keyword will be ignored");

  float cpars[6];
  
  for (int i           = 1; i < parser->ntokens; i++)
  {
    if (parser->token[i].isnumber)
      cpars[i-1]       = (work)(parser->token[i].value);
    else
      Bp3Error("Input::cellkeyword", "A numerical value is expected for CELL");
  }
  
  if (parser->ntokens == 4)
  {
    Bp3Warning("Input::cellkeyword", "Setting alpha = beta = gamma = 90.0");
    cpars[3]           = cpars[4] = cpars[5] = 90.0;
  }
  xtal.cell[nxtal].Setcell(cpars);
}

void Input::centkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 1)
    Bp3Warning("Input::centkeyword", "Extra text in CENTric keyword will be ignored");
  xtal.Setonlycentrics(true);
}

void Input::acenkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens   > 1)
    Bp3Warning("Input::acenkeyword", "Extra text in ACENTric keyword will be ignored");
  xtal.Setonlyacentrics(true);
}

void Input::iampkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens   > 1)
    Bp3Warning("Input::iampkeyword", "Extra text in IAMPlitude keyword will be ignored");
  inputamplitudes       = true;
}

void Input::iintkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens   > 1)
    Bp3Warning("Input::iintkeyword", "Extra text in IINTensity keyword will be ignored");
  inputintensities      = true;
}

void Input::invekeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens   > 1)
    Bp3Warning("Input::invekeyword", "Extra text in INVErt keyword will be ignored");
  xtal.Setinvert(true);
}

void Input::scalkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens   > 1)
    Bp3Warning("Input::scalkeyword", "Extra text in SCALe keyword will be ignored");
  xtal.Setrscale(true);
}

void Input::pluzkeyword()
{
  if (!crystal)
    Bp3Error("Input::pluzkeyword", "XTAL keyword must be given before PLUZ");

  vector<double> userpdluz;
  unsigned pindex(0);
  
  if (parser->token[1].isnumber)
    if (parser->token[1].value         >= 1)
      pindex                            = (unsigned)(parser->token[1].value);
    else
      Bp3Error("Input::pluzkeyword", "A natural number is expected for the index of PLUZ");
  else 	
    Bp3Error("Input::pluzkeyword", "A natural number is expected for the index of PLUZ");

  for (int i                            = 2; i < parser->ntokens-1; i++)
  {
    if (parser->token[i].isnumber)
      userpdluz.push_back((double)(parser->token[i].value));
    else
      Bp3Error("Input::pluzkeyword", "A numerical value is expected for PLUZ");
  }

  // the last token might be the NORE keyword
  
  if (parser->token[parser->ntokens-1].isnumber)
  {
    userpdluz.push_back((double)(parser->token[parser->ntokens-1].value));
    xtal.Setrefinep(pindex-1,true);
  }
  else
    if (!CCP4::ccp4_keymatch("NORE", parser->token[parser->ntokens-1].word))
      Bp3Error("Input::pluzkeyword", "The NOREfine keyword was expected");
    else
      xtal.Setrefinep(pindex-1,false);

  // set or ensure nbins is the same as number of parameters inputted
  if (userpdluz.size()                  > 1)
  {
    if (!xtal.sf[ndata].Getnbins())
      xtal.sf[ndata].Setusernbins(userpdluz.size());
    else if (xtal.sf[ndata].Getnbins() != userpdluz.size())
      Bp3Error("Input::pluzkeyword",
	       "number of BINS does not match PLUZ parameters inputted");
    xtal.Setpluz(pindex-1,userpdluz);

    for (unsigned p                     = 0;   p < pindex; p++)
      for (unsigned r                   = p+1; r < pindex; r++)
	if (xtal.Getpluz(p).size())
	  if (xtal.Getpluz(r).size())
	    if (xtal.Getpluz(p).size() != xtal.Getpluz(r).size())
	      Bp3Error("Input::pluzkeyword",
		       "number of BINS does not match for two different PLUZ inputted");
	      
  }
  else
    xtal.Setdefaultpluz(pindex-1,userpdluz[0]);

  xtal.Setuserpluz(pindex-1,true);
}

// DATA KEYWORDS

bool Input::datakeywords()
{
  if (CCP4::ccp4_keymatch("DNAM", parser->keyword))
    dnamkeyword();
  else if (CCP4::ccp4_keymatch("MTZI", parser->keyword))
    mtzikeyword();
  else if (CCP4::ccp4_keymatch("SCAI", parser->keyword))
    scaikeyword();
  else if (CCP4::ccp4_keymatch("UNKF", parser->keyword))
    unkfkeyword();
  else if (CCP4::ccp4_keymatch("COLU", parser->keyword))
    colukeyword();
  else if (CCP4::ccp4_keymatch("FORM", parser->keyword))
    formkeyword();
  else if (CCP4::ccp4_keymatch("RESO", parser->keyword))
    resokeyword();
  else if (CCP4::ccp4_keymatch("BINS", parser->keyword))
    binskeyword();
  else if (CCP4::ccp4_keymatch("KSCA", parser->keyword))
    kscalekeyword();
  else if (CCP4::ccp4_keymatch("BSCA", parser->keyword))
    bscalekeyword();
  else if (CCP4::ccp4_keymatch("ISOE", parser->keyword))
    isoekeyword();
  else if (CCP4::ccp4_keymatch("ANOE", parser->keyword))
    anoekeyword();
  else if (CCP4::ccp4_keymatch("CORE", parser->keyword))
    corekeyword();
  else if (CCP4::ccp4_keymatch("SDLU", parser->keyword))
    sdlukeyword();
  else if (CCP4::ccp4_keymatch("EXCL", parser->keyword))
    exclkeyword();
  else
    return false;

  return true;
}

void Input::dnamkeyword()
{
  if (nxtal             < 0)
    Bp3Error("Input::dnamkeyword", "No XTAL has been defined");
    
  dataset               = true;
  ndata++;
  nwave++;
  xtal.sf.push_back(Sfdata(mbin,refbin));
  xtal.sf[ndata].Setnxtal(nxtal);
  xtal.sf[ndata].Setnwave(nwave);
  
  // check to see if we have no more than two tokens
  if (parser->ntokens   > 2)
    Bp3Warning("Input::dnamkeyword", "Extra text in DNAMe keyword will be ignored");

  if (parser->token[1].isstring)
    xtal.sf[ndata].Setname(parser->token[1].fullstring);
  else
    Bp3Error("Input::dnamkeyword", "A string is expected for the DNAMe's name");
    
}

void Input::mtzikeyword()
{
  if (!dataset)
    Bp3Error("Input::mtzikeyword","DNAMe keyword must be given before MTZIn");  

  // check to see if we have no more than two tokens
  if (parser->ntokens   > 2)
    Bp3Warning("Input::mtzikeyword", "Extra text in MTZIn keyword will be ignored");

  if (parser->token[1].isstring)
    xtal.sf[ndata].Setmtzin(parser->token[1].fullstring);
  else
    Bp3Error("Input::dnamkeyword", "A string is expected after MTZIn");
}

void Input::scaikeyword()
{
  if (!dataset)
    Bp3Error("Input::scaikeyword","DNAMe keyword must be given before SCAIn");  

  // check to see if we have no more than two tokens
  if (parser->ntokens   > 2)
    Bp3Warning("Input::scaikeyword", "Extra text in SCAIn keyword will be ignored");

  if (parser->token[1].isstring)
    xtal.sf[ndata].Setscain(parser->token[1].fullstring);
  else
    Bp3Error("Input::scaikeyword", "A string is expected after SCAIn");
}

void Input::unkfkeyword()
{
  if (!dataset)
    Bp3Error("Input::unkfkeyword","DNAMe keyword must be given before UNKFile");  

  // check to see if we have no more than two tokens
  if (parser->ntokens   > 2)
    Bp3Warning("Input::unkfkeyword", "Extra text in UNKFile keyword will be ignored");

  if (parser->token[1].isstring)
    unknownfiletype.push_back(parser->token[1].fullstring);
  else
    Bp3Error("Input::unkfkeyword", "A string is expected after UNKFile");
}

void Input::colukeyword()
{
  if (!dataset)
    Bp3Error("Input::colukeyword","DNAMe keyword must be given before COLUmn");  
  
  bool gotd(false);
  bool gotsig(false);
  bool gotdp(false);
  bool gotsigp(false);
  bool gotdm(false);
  bool gotsigm(false);
  
  for (int i                 = 1; i < parser->ntokens; i++)
  {
    if ((CCP4::ccp4_keymatch("F",parser->token[i].fullstring) ||
 	 CCP4::ccp4_keymatch("I", parser->token[i].fullstring))
	 && !gotd)
    {
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsdata(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in F or I");
      gotd                   = true;
      if (CCP4::ccp4_keymatch("F",parser->token[i].fullstring))
	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "AMPLITUDE") )
	  xtal.sf[ndata].Settype("AMPLITUDE");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
      else
	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "INTENSITY") )
	  xtal.sf[ndata].Settype("INTENSITY");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
    }
    else if ((CCP4::ccp4_keymatch("SF", parser->token[i].fullstring) ||
	      CCP4::ccp4_keymatch("SI", parser->token[i].fullstring)) && !gotsig)
    {
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsdev(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in SF or SI");
      gotsig                 = true;
      if (CCP4::ccp4_keymatch("SF",parser->token[i].fullstring))
 	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "AMPLITUDE") )
 	  xtal.sf[ndata].Settype("AMPLITUDE");
 	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
      else
	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "INTENSITY") )
	  xtal.sf[ndata].Settype("INTENSITY");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
    }
    if ((CCP4::ccp4_keymatch("F+",parser->token[i].fullstring) ||
	 CCP4::ccp4_keymatch("I+",parser->token[i].fullstring)) && !gotdp)
    {
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsdatap(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in F+ or I+");
      if (CCP4::ccp4_keymatch("F+",parser->token[i].fullstring))
 	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "AMPLITUDE") )
	  xtal.sf[ndata].Settype("AMPLITUDE");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
      else
	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "INTENSITY") )
	  xtal.sf[ndata].Settype("INTENSITY");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
      gotdp                  = true;
    }
    else if ((CCP4::ccp4_keymatch("SF+", parser->token[i].fullstring) ||
	      CCP4::ccp4_keymatch("SI+", parser->token[i].fullstring)) && !gotsigp)
    {
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsdevp(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in SF+ or SI+");
      gotsigp                = true;
      if (CCP4::ccp4_keymatch("SF+",parser->token[i].fullstring))
 	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "AMPLITUDE") )
	  xtal.sf[ndata].Settype("AMPLITUDE");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
      else
	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "INTENSITY") )
	  xtal.sf[ndata].Settype("INTENSITY");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
    }
    else if ((CCP4::ccp4_keymatch("F-", parser->token[i].fullstring) ||
	      CCP4::ccp4_keymatch("I-", parser->token[i].fullstring)) && !gotdm )
    {
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsdatam(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in F- or I-");
      gotdm                  = true;
      if (CCP4::ccp4_keymatch("F-",parser->token[i].fullstring))
 	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "AMPLITUDE") )
	  xtal.sf[ndata].Settype("AMPLITUDE");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
      else
	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "INTENSITY") )
	  xtal.sf[ndata].Settype("INTENSITY");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
    }
    else if ((CCP4::ccp4_keymatch("SF-", parser->token[i].fullstring) ||
 	      CCP4::ccp4_keymatch("SI-", parser->token[i].fullstring)) && !gotsigm )
    {
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsdevm(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in SF- or SI-");
      gotsigm                = true;
      if (CCP4::ccp4_keymatch("SF-",parser->token[i].fullstring))
 	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "AMPLITUDE") )
	  xtal.sf[ndata].Settype("AMPLITUDE");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
      else
	if ( (xtal.sf[ndata].Gettype() == "") ||
	     (xtal.sf[ndata].Gettype() == "INTENSITY") )
	  xtal.sf[ndata].Settype("INTENSITY");
	else
	  Bp3Error("Input::colukeyword","Set either amplitudes or intensities");
    }
    else if (CCP4::ccp4_keymatch("PHIB",parser->token[i].fullstring))
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsphi(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in PHIB");
    else if (CCP4::ccp4_keymatch("FOM",parser->token[i].fullstring))
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setsfom(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in FOM");
    else if (CCP4::ccp4_keymatch("HLA",parser->token[i].fullstring))
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setshla(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in HLA");
    else if (CCP4::ccp4_keymatch("HLB",parser->token[i].fullstring))
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setshlb(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in HLB");
    else if (CCP4::ccp4_keymatch("HLC",parser->token[i].fullstring))
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setshlc(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in HLC");
    else if (CCP4::ccp4_keymatch("HLD",parser->token[i].fullstring))
      if ((i+1)              < parser->ntokens)
	xtal.sf[ndata].Setshld(parser->token[i+1].fullstring);
      else
	Bp3Error("Input::colukeyword","Parsing error in HLD");
    else
      for (unsigned m          = 1; m < 100; m++)
      {
        char number[4];
        sprintf(number,"%u",m);
        string fc              = number;
	string pc              = "PC" + fc;
        fc                     = "FC" + fc;	
	if (CCP4::ccp4_keymatch(fc.c_str(),parser->token[i].fullstring))
	  if ((i+1)            < parser->ntokens)
	    xtal.sf[ndata].Setsfmodel(m,parser->token[i+1].fullstring);
	  else
	    Bp3Error("Input::colukeyword","Parsing error in FC");
	else if (CCP4::ccp4_keymatch(pc.c_str(),parser->token[i].fullstring))
	  if ((i+1)            < parser->ntokens)
	    xtal.sf[ndata].Setspmodel(m,parser->token[i+1].fullstring);
	  else
	    Bp3Error("Input::colukeyword","Parsing error in PHIC");
      }
  }
}

void Input::formkeyword()
{
  if (!dataset)
    Bp3Error("Input::formkeyword","DNAMe keyword must be given before FORM");

  if (parser->ntokens > 6)
    Bp3Warning("Input::formkeyword",
	       "Extra text in FORM keyword will be ignored");

  if (parser->ntokens < 6)
    Bp3Error("Input::formkeyword",
	     "Please specify both FP and FPP in FORM keyword");

  // atomindex
  unsigned a(0);
  bool formset            = false;
  string localatomname;
  
  if (parser->token[1].isstring)
    localatomname         = parser->token[1].fullstring;
  else
    Bp3Error("Input::formkeyword", "A string is expected for the atom name");

  // check to see if the atom's form factors have already been stored

  unsigned nx             = (unsigned) nxtal;
  unsigned nw             = (unsigned) nwave;

  for (unsigned i         = 0; i < localatomname.size(); i++)
    localatomname[i]      = toupper(localatomname[i]);
  
  for (unsigned f         = 0; f < mdl.form.size(); f++)
    if ((localatomname    == mdl.form[f].Getname()) && (nx == mdl.form[f].Getcrystal()))
    {
      formset             = true;
      a                   = f;
    }
  
  if (!formset)
  {
    mdl.form.push_back(Scatter(localatomname, nx));
    a                     = mdl.form.size() - 1;
  }

  mdl.form[a].Setnwave(nw+1);
  
  for (int i              = 1; i < parser->ntokens; i++)
  {
    if (CCP4::ccp4_keymatch("FP", parser->token[i].word))
      if ((i+1)           < parser->ntokens)
      {
	if (parser->token[i+1].isnumber)
	  mdl.form[a].Setfp(nw, (double) parser->token[i+1].value);
	else
	  Bp3Error("Input::formkeyword", "A numerical value is expected for FP");
      }
      else
	Bp3Error("Input::formkeyword", "Parsing error");

    if (CCP4::ccp4_keymatch("FPP", parser->token[i].word))
      if ((i+1)           < parser->ntokens)
      {
	if (parser->token[i+1].isnumber)
	  mdl.form[a].Setfpp(nw, (double) parser->token[i+1].value);
	else
	  Bp3Error("Input::formkeyword", "A numerical value is expected for FPP");
      }
      else
	Bp3Error("Input::formkeyword", "Parsing error"); 
  }
}

void Input::resokeyword()
{
  if (!dataset)
    Bp3Error("Input::resokeyword","DNAMe keyword must be given before RESOlution");

  work limit[2];

  if (parser->token[1].isnumber)
  {
    for (unsigned i     = 1; i < 3; i++)
      if (parser->token[i].isnumber)
	limit[i-1]      = (work) parser->token[i].value;
      else
	Bp3Error("Input::resokeyword", "A numerical value is expected for the resolution limit - and both a high and low resolution limit must be given");
	
    xtal.sf[ndata].Sethires(std::min(limit[0],limit[1]));
    xtal.sf[ndata].Setlowres(std::max(limit[0],limit[1]));

    if (parser->ntokens > 3)
      Bp3Error("Input::resokeyword","Parsing error");
  }
  else if (parser->token[1].isstring)
    Bp3Error("Input::resokeyword","Parsing error");
  else
    Bp3Error("Input::resokeyword","RESOlution limits is expected after the RESO keyword");
}

void Input::binskeyword()
{
  if (!dataset)
    Bp3Error("Input::binskeyword","DNAMe keyword must be given before BINS");

  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Input::binskeyword", "Extra text in BINS keyword will be ignored");

  if (!parser->token[1].isnumber)	
    Bp3Error("Input::binskeyword", "A numerical value is expected for BINS");
  
  xtal.sf[ndata].Setusernbins((unsigned)(parser->token[1].value));
}

void Input::kscalekeyword()
{
  if (!dataset)
    Bp3Error("Input::kscalekeyword","DNAMe keyword must be given before KSCAle");

  if (!parser->token[1].isnumber)	
    Bp3Error("Input::kscalekeyword", "A numerical value is expected for KSCAle");

  xtal.sf[ndata].Setkscale((double)(parser->token[1].value));

  // check to see if we have no more than one token
  if (parser->ntokens  > 3)
    Bp3Warning("Input::kscalekeyword",
	       "Extra text in KSCAle keyword will be ignored");

  if (parser->ntokens == 3)
    if (!CCP4::ccp4_keymatch("REFI", parser->token[2].word))
      Bp3Error("Input::kscalekeyword", "The REFIne keyword was expected");
    else
      xtal.sf[ndata].Setrefinek(true);
}

void Input::bscalekeyword()
{
  if (!dataset)
    Bp3Error("Input::bscalekeyword","DNAMe keyword must be given before BSCAle");

  if (!parser->token[1].isnumber)	
    Bp3Error("Input::bscalekeyword", "A numerical value is expected for BSCAle");

  xtal.sf[ndata].Setbiso((double)(parser->token[1].value));

  // check to see if we have no more than one token
  if (parser->ntokens > 3)
    Bp3Warning("Input::bscalekeyword",
	       "Extra text in BSCAle keyword will be ignored");

  if (!CCP4::ccp4_keymatch("REFI", parser->token[2].word))
    Bp3Error("Input::bscalekeyword", "The REFIne keyword was expected");
  else
    xtal.sf[ndata].Setrefineb(true);
}

void Input::isoekeyword()
{
  if (!dataset)
    Bp3Error("Input::isoekeyword", "DNAMe keyword must be given before ISOE");

  vector<double> userdluz;

  if (parser->token[1].isstring)
  {
    if (!CCP4::ccp4_keymatch("NORE", parser->token[parser->ntokens-1].word))
      Bp3Error("Input::isoekeyword", "The NORE keyword was expected");
    else
      xtal.sf[ndata].Setrefined(false);
  }
  else
  {
    for (int i = 1; i < parser->ntokens-1; i++)
    {
      if (parser->token[i].isnumber)
	userdluz.push_back((double)(parser->token[i].value));
      else
	Bp3Error("Input::isoekeyword", "A numerical value is expected for ISOE");
    }

    // the last token might be the NORE keyword
    if (parser->token[parser->ntokens-1].isnumber)
      userdluz.push_back((double)(parser->token[parser->ntokens-1].value));
     else
       if (!CCP4::ccp4_keymatch("NORE", parser->token[parser->ntokens-1].word))
	 Bp3Error("Input::isoekeyword", "The NORE keyword was expected");
       else
	 xtal.sf[ndata].Setrefined(false);

    // set or ensure nbins is the same as number of parameters inputted

    if (!xtal.sf[ndata].Getnbins())
      xtal.sf[ndata].Setusernbins(userdluz.size());
    else if (xtal.sf[ndata].Getnbins() != userdluz.size())
      Bp3Error("Input::isoekeyword",
	       "number of BINS does not match ISOE parameters inputted");

    xtal.sf[ndata].Setdluz(userdluz);
  }
}

void Input::anoekeyword()
{
  if (!dataset)
    Bp3Error("Input::anoekeyword", "DNAMe keyword must be given before ANOE");

  vector<double> useradluz;
  
  for (int i = 1; i < parser->ntokens-1; i++)
  {
    if (parser->token[i].isnumber)
      useradluz.push_back((double)(parser->token[i].value));
    else
      Bp3Error("Input::anoekeyword", "A numerical value is expected for ANOE");
  }

  // the last token might be the NORE keyword
  
  if (parser->token[parser->ntokens-1].isnumber)
    useradluz.push_back((double)(parser->token[parser->ntokens-1].value));
  else
    if (!CCP4::ccp4_keymatch("NORE", parser->token[parser->ntokens-1].word))
      Bp3Error("Input::anoekeyword", "The NOREFine keyword was expected");
    else
      xtal.sf[ndata].Setrefinead(false);

  // set or ensure nbins is the same as number of parameters inputted
  if (!xtal.sf[ndata].Getnbins())
    xtal.sf[ndata].Setusernbins(useradluz.size());
  else if (xtal.sf[ndata].Getnbins() != useradluz.size())
    Bp3Error("Input::anoekeyword",
	     "number of BINS does not match ANOE parameters inputted");
      
  xtal.sf[ndata].Setadluz(useradluz);
}

void Input::corekeyword()
{
  if (!dataset)
    Bp3Error("Input::corekeyword", "DNAMe keyword must be given before CORE");

  vector<double> usereluz;
  
  for (int i = 1; i < parser->ntokens-1; i++)
  {
    if (parser->token[i].isnumber)
      usereluz.push_back((double)(parser->token[i].value));
    else
      Bp3Error("Input::corekeyword", "A numerical value is expected for CORE");
  }

  // the last token might be the NORE keyword
  
  if (parser->token[parser->ntokens-1].isnumber)
    usereluz.push_back((double)(parser->token[parser->ntokens-1].value));
  else
    if (!CCP4::ccp4_keymatch("REFI", parser->token[parser->ntokens-1].word))
      Bp3Error("Input::corekeyword", "The NOREfine keyword was expected");
    else
      xtal.sf[ndata].Setrefinee(true);

  // set or ensure nbins is the same as number of parameters inputted
  if (!xtal.sf[ndata].Getnbins())
    xtal.sf[ndata].Setusernbins(usereluz.size());
  else if (xtal.sf[ndata].Getnbins() != usereluz.size())
    Bp3Error("Input::corekeyword",
	     "number of BINS does not match CORE parameters inputted");
      
  xtal.sf[ndata].Seteluz(usereluz);
}

void Input::sdlukeyword()
{
  if (!dataset)
    Bp3Error("Input::sdlukeyword", "DNAMe keyword must be given before SDLUz");

  vector<double> usersdluz;
  
  for (int i = 1; i < parser->ntokens-1; i++)
  {
    if (parser->token[i].isnumber)
      usersdluz.push_back((double)(parser->token[i].value));
    else
      Bp3Error("Input::sdlukeyword", "A numerical value is expected for SDLUz");
  }

  // the last token might be the NORE keyword
  
  if (parser->token[parser->ntokens-1].isnumber)
    usersdluz.push_back((double)(parser->token[parser->ntokens-1].value));
  else
    if (!CCP4::ccp4_keymatch("NORE", parser->token[parser->ntokens-1].word))
      Bp3Error("Input::sdlukeyword", "The NOREfine keyword was expected");
    else
      xtal.sf[ndata].Setrefinesd(false);

  // set or ensure nbins is the same as number of parameters inputted
  if (usersdluz.size() > 1)
  {
    if (!xtal.sf[ndata].Getnbins())
      xtal.sf[ndata].Setusernbins(usersdluz.size());
    else if (xtal.sf[ndata].Getnbins() != usersdluz.size())
      Bp3Error("Input::sdlukeyword",
	       "number of BINS does not match SDLUz parameters inputted");
      
    xtal.sf[ndata].Setsdluz(usersdluz);
  }
  else
    xtal.sf[ndata].Setdefaultsdluz(usersdluz[0]);
}

void Input::exclkeyword()
{
  if (!dataset)
    Bp3Error("Input::exclkeyword","DNAMe keyword must be given before EXCLude");  
  
  for (int i                 = 1; i < parser->ntokens; i++)
  {
    if (CCP4::ccp4_keymatch("SIGF",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	if (parser->token[i+1].isnumber)
	  xtal.sf[ndata].Setsigmacut(parser->token[i+1].value);
	else
	  Bp3Error("Input::exclkeyword","A number is expected for the sigma cut-off");
      else
	Bp3Error("Input::exclkeyword","Parsing error");
    }
    else if (CCP4::ccp4_keymatch("SANO",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	if (parser->token[i+1].isnumber)
	  xtal.sf[ndata].Setanosigmacut(parser->token[i+1].value);
	else
	  Bp3Error("Input::exclkeyword","A number is expected for the ano sigma cut-off");
      else
	Bp3Error("Input::exclkeyword","Parsing error");
    } 
    else if (CCP4::ccp4_keymatch("SISO",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	if (parser->token[i+1].isnumber)
	  xtal.sf[ndata].Setisosigmacut(parser->token[i+1].value);
	else
	  Bp3Error("Input::exclkeyword","A number is expected for the iso sigma cut-off");
      else
	Bp3Error("Input::exclkeyword","Parsing error");
    }
  }
}

// MODEL KEYWORDS

bool Input::modelkeywords()
{
  if (CCP4::ccp4_keymatch("MODL", parser->keyword))       // model keywords
    modlkeyword();
  else if (CCP4::ccp4_keymatch("SITE", parser->keyword))            
    sitekeyword();
  else if (CCP4::ccp4_keymatch("ATOM", parser->keyword))       
    atomkeyword();
  else if (CCP4::ccp4_keymatch("XYZ", parser->keyword))
    xyzkeyword();
  else if (CCP4::ccp4_keymatch("NUMB", parser->keyword))
    numbkeyword();  
  else if (CCP4::ccp4_keymatch("OCCU", parser->keyword))
    occukeyword();
  else if (CCP4::ccp4_keymatch("BISO", parser->keyword))
    bisokeyword();
  else if (CCP4::ccp4_keymatch("UANO", parser->keyword))
    uanokeyword();
  else if (CCP4::ccp4_keymatch("NORE", parser->keyword))
    norekeyword();
  else
    return false;
  
  return true;
}

void Input::sitekeyword()
{
  if (setwithmodl || setwithatom || setnumb)
    Bp3Error("Input::sitekeyword", "Do not use SITE keyword with MODL, XYZ, or NUMB keyword");

  if (parser->ntokens       < 5)
    Bp3Error("Input::sitekeyword", "Site number and X, Y and Z values are needed");
  else if (parser->ntokens >= 5)
  {    
    unsigned sitenumber(0);

    if (parser->token[1].isnumber)
      sitenumber            = (unsigned)parser->token[1].value;
    else
      Bp3Error("Input::sitekeyword", "A numerical value is expected for the SITE number");
      
    vector<double> fracx(3, ZERO);
    bool refx(true), refy(true), refz(true);
    
    // get the atomic coordinates
    for (unsigned i         = 2; i < 5; i++) 
    {
      if (parser->token[i].isnumber)
	fracx[i-2]          = (double)parser->token[i].value;
      else
	Bp3Error("Input::sitekeyword", "A numerical value is expected for the fractional coordinates");
    }
    
    // check if user does not wish to refine x, y and/or z
    if (parser->ntokens     > 5)
      if (!CCP4::ccp4_keymatch("NORE", parser->token[5].word))
	Bp3Error("Input::sitekeyword", "The NORE keyword was expected");
      for (int i            = 6; i < parser->ntokens; i++)
      {
	if (parser->token[i].word[0] == 'X')
	  refx              = false;
	if (parser->token[i].word[0] == 'Y')
	  refy              = false;
	if (parser->token[i].word[0] == 'Z')
	  refz              = false;
      }
      
      if (mdl.site.size()  == (sitenumber-1))
	mdl.site.push_back(Site(mdl.site.size(), fracx[0], fracx[1], fracx[2], refx, refy, refz));
      else
	Bp3Error("Input::sitekeyword", "Incorrect site number");

      setwithsite           = true;
  }
}

void Input::atomkeyword()
{
  if (atomready())
    atomset();

  // increment atom counter
  natom++;
  // reset atomic parameter
  resetatomin();
  
  if (parser->ntokens   < 2)
    Bp3Error("Input::atomkeyword", "A name is expected after the ATOM keyword");

  if (parser->token[1].isstring)
    atomname            = parser->token[1].fullstring;
  else
    Bp3Error("Input::atomkeyword", "A string is expected for ATOM");

  for (unsigned i       = 0; i < atomname.size(); i++)
    atomname[i]         = toupper(atomname[i]);
  
  if (parser->ntokens   > 2)
    if (setwithsite)
      if (CCP4::ccp4_keymatch("SITE", parser->token[2].fullstring))
      {
	if (parser->token[3].isnumber)
	  nsite         = (unsigned)(parser->token[3].value - 1);
	else
	  Bp3Error("Input::atomkeyword", "A numerical value is expected for SITE");

	setsite         = true;
      }
      else
	Bp3Error("Input::atomkeyword", "SITE subkeyword was expected");
    else
      Bp3Error("Input::atomkeyword", "SITE has not been set");

}

void Input::modlkeyword()
{
  if (setwithsite || setwithatom || setnumb)
    Bp3Error("Input::xyzkeyword", "Do not use MODL keyword with SITE, XYZ or MODL keyword");

  bfacset             = false;
  if (parser->ntokens > 4)
    Bp3Warning("Input:modlkeyword", "Extra text after MODL keyword will be ignored");

  if (nxtal           < 0)
    Bp3Error("Input::modlkeyword", "The XTAL keyword must be given before MODL");

  if (parser->ntokens > 2)
    for (int i        = 1; i < parser->ntokens; i++)
      if (CCP4::ccp4_keymatch("BFAC", parser->token[i].word))
	if (parser->token[i+1].isnumber)
	{
	  bfac.resize(1,parser->token[i+1].value);
	  bfacset     = true;
	}
    
  if (parser->token[1].isstring)
    readpdb(parser->token[1].fullstring);
  else
    Bp3Error("Input::modlkeyword", "A string is expected for the MODL keyword");
  setwithmodl         = true;
}

void Input::xyzkeyword()
{
  if (setwithsite || setwithmodl || setnumb)
    Bp3Error("Input::xyzkeyword", "Do not use XYZ keyword with SITE, MODL or NUMB keyword");

  if (atomname             == "")
    Bp3Error("Input::xyzkeyword", "Must give ATOM keyword before XYZ");

  nsite++;
  
  vector<double> fracx(3, ZERO);

  bool refx(true), refy(true), refz(true);

  if (parser->ntokens       < 4)
    Bp3Error("Input::xyzkeyword", "X, Y and Z values are needed");
  else if (parser->ntokens >= 4)
  {    
    // get the atomic coordinates
    for (unsigned i         = 1; i < 4; i++)
    {
      if (parser->token[i].isnumber)
	fracx[i-1]          = (double)parser->token[i].value;
      else
	Bp3Error("Input::xyzkeyword", "A numerical value is expected for the fractional coordinates");
    }
    
    // check if user does not wish to refine x, y and/or z
    if (parser->ntokens     > 4) 
    {
      if (!CCP4::ccp4_keymatch("NORE", parser->token[4].word))
	Bp3Error("Input::xyzkeyword", "The NOREfine keyword was expected");
      for (int i            = 5; i < parser->ntokens; i++)
      {
	if ( (parser->token[i].word[0] == 'X') || (parser->token[i].word[0] == 'x') )
	  refx              = false;
	if ( (parser->token[i].word[0] == 'Y') || (parser->token[i].word[0] == 'y') )
	  refy              = false;
	if ( (parser->token[i].word[0] == 'Z') || (parser->token[i].word[0] == 'z') )
	  refz              = false;
      }
    }
  }
  
  // set the site
  mdl.site.push_back(Site(nsite, fracx[0], fracx[1], fracx[2], refx, refy, refz));
  setwithatom               = true;
  setsite                   = true;
}

void Input::numbkeyword()
{
  if (setwithsite || setwithmodl || setwithatom)
    Bp3Error("Input::numbkeyword", "Do not use NUMB keyword with SITE, MODL or XYZ keyword");

  if (atomname       == "")
    Bp3Error("Input::numbkeyword", "Must give ATOM keyword before NUMB");

  /*
  if (mdl.atom.size() > 0)
    Bp3Error("Input::numbkeyword", "Can not set NUMB and also assign xyz coordinates with MODL or XYZ");
  */
  
  nsite++;  

  if (parser->ntokens > 2)
    if (CCP4::ccp4_keymatch("NORE", parser->token[2].word))
      refnumb         = false;
    else
      Bp3Error("Input::occukeyword", "The NOREfine keyword was expected");

  // let occupancy hold expected number of atoms
  if (parser->token[1].isnumber)
  {
    number            = (double)parser->token[1].value;
    occ               = ONE;
  }
  else
    Bp3Error("Input::numbkeyword", "A numerical value is expected for the number of expected atoms");
  
  // ***NSP - don't need to set the site
  // mdl.site.push_back(Site(nsite, 0, fracx[1], fracx[2], refx, refy, refz));
 
  setnumb             = true; 
  setsite             = seto = setb = true;

  /*
  if (atomready())
    atomset();
  */
}

void Input::occukeyword()
{
  if (atomname           == "")
    Bp3Error("Input::occukeyword", "Must give ATOM keyword before OCCU");

  if (setnumb)
    Bp3Error("Input::occukeyword", "Can not set NUMB keyword with OCCU keyword");

  // get the occupancy
  if (parser->token[1].isnumber)
  {
    occ               = (double)parser->token[1].value;
    seto              = true;
  }
  else
    Bp3Error("Input::occukeyword", "A numerical value is expected for OCCU");

  // check if user does not wish to refine
  if (parser->ntokens > 2)
    if (CCP4::ccp4_keymatch("NORE", parser->token[2].word))
      refo            = false;
    else
      Bp3Error("Input::occukeyword", "The NOREfine keyword was expected");

  /*
  if (atomready())
    atomset();
  */
}

void Input::bisokeyword()
{
  if (atomname           == "")
    Bp3Error("Input::bisokeyword", "Must give ATOM keyword before BISO");

  // get the bfactor
  if (parser->token[1].isnumber)
  {
    bfac.resize(1);
    bfac[0]           = (double)parser->token[1].value;
    setb              = true;
    biso              = true;
  }
  else
    Bp3Error("Input::bisokeyword", "A numerical value is expected for BISO");

  // check if user does not wish to refine
  if (parser->ntokens > 2)
    if (CCP4::ccp4_keymatch("NORE", parser->token[2].word))
      refb            = false;
    else
      Bp3Error("Input::bisokeyword", "The NOREfine keyword was expected");

  if (atomready())
    atomset();
}

void Input::uanokeyword()
{

  if (atomname           == "")
    Bp3Error("Input::uanokeyword", "Must give ATOM keyword before UANO");
    
  // get the anisotropic bfactor
    biso                  = false;
  if (parser->token[1].isnumber)
  {
    bfac[0]               = (double)parser->token[1].value;
    setb                  = true;
  }
  else
    Bp3Error("Input::uanokeyword", "A numerical value is expected for the BFACtor");

  if (parser->ntokens     > 2)
  {
    if (parser->token[2].isnumber)
    {
      if (parser->ntokens != 7)
	Bp3Error("Input::uanokeyword", "1 OR 6 numerical values are needed for BANO");
      bfac.resize(6);
      for (unsigned i     = 2; i < 7; i++)
	if (parser->token[i].isnumber)
	  bfac[i-1]       = (double)parser->token[i].value;
	else
	  Bp3Error("Input::uanokeyword", "A numerical value is expected for the bfactor");
      if (parser->ntokens > 7)
	if (CCP4::ccp4_keymatch("NORE", parser->token[7].word))
	  refb            = false;
	else
	  Bp3Warning("Input::uanokeyword", "Extra text in BANO keyword will be ignored");      
    }
    else if (CCP4::ccp4_keymatch("NORE", parser->token[2].word))
      refb                = false;
    else
      Bp3Error("Input::uanokeyword", "The NORE keyword was expected");
  }
  else
  {
    bfac.resize(1);
    Bp3Warning("Input::uanokeyword", "The isotropic B will be converted to anisotropic");
  }
  
  if (atomready())
    atomset();
}

void Input::norekeyword()
{
  if (parser->ntokens > 1)
  {
    for (int i        = 1; i < parser->ntokens; i++)
      if (CCP4::ccp4_keymatch("XYZ",parser->token[i].word))
	norefxyz      = true;
      else if (CCP4::ccp4_keymatch("OCCU",parser->token[i].word))
 	norefocc      = true;
      else if (CCP4::ccp4_keymatch("BFAC",parser->token[i].word))
 	norefbfac     = true;
      else
	Bp3Error("Input::norekeyword", "Subkeyword of NOREf not recognized");
  }
  else
    Bp3Error("Input::norekeyword", "XYZ, OCCU, and/or BFAC is expected after NOREf");

}

// UTILITY FUNCTION FOR MODEL

void Input::atomset()
{
  bool gotform(false);
  unsigned nform(0);
  
  // check to see if a form factor has already been set for this atom type

  unsigned nx       = (unsigned)nxtal;

  for (unsigned i   = 0; i < mdl.form.size(); i++)
    if ((atomname  == mdl.form[i].Getname()) && (nx == mdl.form[i].Getcrystal()) )
    {
      gotform       = true;
      nform         = i;
    }
  
  // if no form exists for the atom, then create one
  if (!gotform)
  {
    nform           = mdl.form.size();
    mdl.form.push_back(Scatter(atomname,nx));
  }

  unsigned ns       = (unsigned)nsite;

  mdl.atom.push_back(Atom(atomname, ns, nform, nx, biso, bfac, occ,  refb, refo));


  if (setnumb)
  {
    mdl.atom[mdl.atom.size()-1].Setnumber(number);
    if (!refnumb)
      mdl.atom[mdl.atom.size()-1].Setrefinenumb(false);
  }
  
  resetatomin();
}

void Input::readpdb(const char* pdbname)
{
  FILE *pdbfile(fopen(pdbname,"r"));

  if (pdbfile                           == NULL)
  {
    string temp(pdbname);
    Bp3Error("Likelihood::readpdb","File " + temp + " could not be opened for reading");
  }
  else
  {
    xtal.pdbfilename[nxtal]              = pdbname;
    printf("\nCoordinates from pdb read in:\n\n");
    printf(" ATOM   X       Y       Z     OCC    BFAC\n");
    while (!feof(pdbfile))
    {
      const int linesize(90);
      char line[linesize+1]              = "";
      if (fgets(line,linesize,pdbfile)  != NULL)
      {
	int number(-1);
	char atom[5], residue[6], label[7];
	float x, y, z, oc, bfa;
	sscanf(line,"%6s%5d%4s%5s%*c%*c%*c%*c%*c%*c%*c%*c%*c%*c%8f%8f%8f%6f%6f",label, &number, atom, residue, &x, &y, &z, &oc, &bfa);
	if ( ( (strncmp(label,"ATOM",4) == 0) || (strncmp(label,"HETATM",6) == 0) ) && (number != -1)  )
	{
	  string cline                   = line;
	  sscanf(cline.substr(30,8).c_str(),"%8f",&x);
	  sscanf(cline.substr(38,8).c_str(),"%8f",&y);
	  sscanf(cline.substr(46,8).c_str(),"%8f",&z);
	  sscanf(cline.substr(54,6).c_str(),"%6f",&oc);
	  sscanf(cline.substr(60,6).c_str(),"%6f",&bfa);
	  
	  bfac.resize(1);

	  if (!bfacset)
	    bfac[0]                     = (double)bfa;
	  nsite++;
	  printf("%3s%8.3f%8.3f%8.3f%6.2f %6.2f\n",atom,x,y,z,oc,bfa);
	  double dx((double)x), dy((double)y), dz((double)z);
	  occ                           = (double)oc;
	  mdl.site.push_back(Site(mdl.site.size(),dx,dy,dz,true,true,true));
	  atomname                      = atom;
	  for (unsigned i               = 0; i < atomname.size(); i++)
	    if (isalpha(atomname[i]))
		atomname[i]             = toupper(atomname[i]);
	  setsite                       = biso = refb = refo = true;
	  frac[nxtal]                   = false;
	  atomset();
	}
      }
    }
    fclose(pdbfile);

    printf("\n");
    bfacset                   = false;
  }
}


// MINIMIZER KEYWORDS

bool Input::minimizerkeywords()
{
  if (CCP4::ccp4_keymatch("CYCL", parser->keyword))
    cyclkeyword();
  else if (CCP4::ccp4_keymatch("SEAR", parser->keyword))            
    searkeyword();
  else if (CCP4::ccp4_keymatch("LINE", parser->keyword))            
    linekeyword();
  else if (CCP4::ccp4_keymatch("NORM", parser->keyword))            
    normkeyword();
  else if (CCP4::ccp4_keymatch("TPAR", parser->keyword))            
    tparkeyword();
  else if (CCP4::ccp4_keymatch("WOLF", parser->keyword))
    wolfkeyword();
  else
    return false;

  return true;
}

void Input::cyclkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::cyclkeyword", "Extra text in CYCLe keyword will be ignored");

  if (parser->token[1].isnumber)
    mini.Setiterations((unsigned)parser->token[1].value);
  else
    Bp3Error("Input::cyclkeyword", 
	     "A numerical value is expected for the CYCle number");
}

void Input::searkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::searkeyword", "Extra text in SEARch keyword will be ignored");

  if (parser->token[1].isstring)
    mini.Setsearch(parser->token[1].fullstring);
  else
    Bp3Error("Input::searkeyword", "A string is expected for the SEARch keyword");
}

void Input::linekeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::linekeyword", "Extra text in LINE keyword will be ignored");

  if (parser->token[1].isstring)
    mini.Setline(parser->token[1].fullstring);
  else
    Bp3Error("Input::linekeyword", "A string is expected for the LINE keyword");
}

void Input::normkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::normkeyword", "Extra text in NORM keyword will be ignored");
  
  if (parser->token[1].isnumber)
    mini.Setnormtol(parser->token[1].value);
  else
    Bp3Error("Input::normkeyword", "A number is expected for the NORM keyword");
}

void Input::tparkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::tparkeyword", "Extra text in TPAR keyword will be ignored");
  
  if (parser->token[1].isnumber)
    mini.Setpartol(parser->token[1].value);
  else
    Bp3Error("Input::tparkeyword", "A number is expected for the TPAR keyword");
}

void Input::wolfkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 3)
    Bp3Warning("Input::wolfkeyword", "Extra text in WOLF keyword will be ignored");

  if (parser->ntokens  < 3)
    Bp3Error("Input::wolfkeyword", "alpha and beta values are expected after wolfe");

  if (parser->token[1].isnumber)
    mini.Setalpha(parser->token[1].value);
  else
    Bp3Error("Input::wolfkeyword", "A number is expected for alpha");

  if (parser->token[2].isnumber)
    mini.Setbeta(parser->token[2].value);
  else
    Bp3Error("Input::wolfkeyword", "A number is expected for beta");

}

// LIKELIHOOD KEYWORDS

bool Input::likelihoodkeywords()
{
  if (CCP4::ccp4_keymatch("TARG", parser->keyword))     
    targkeyword();
  else if (CCP4::ccp4_keymatch("FILL", parser->keyword)) 
    fillkeyword();
  else if (CCP4::ccp4_keymatch("NODE", parser->keyword)) 
    nodekeyword();
  else if (CCP4::ccp4_keymatch("OUTI", parser->keyword))            
    outikeyword();
  else if (CCP4::ccp4_keymatch("OUTH", parser->keyword))            
    outhkeyword();
  else if (CCP4::ccp4_keymatch("OUTS", parser->keyword))            
    outskeyword();
  else if (CCP4::ccp4_keymatch("INTE", parser->keyword))            
    intekeyword();
  else if (CCP4::ccp4_keymatch("VERB", parser->keyword))     
    verbkeyword();
  else if (CCP4::ccp4_keymatch("WARN", parser->keyword))     
    warnkeyword();
  else if (CCP4::ccp4_keymatch("TITL", parser->keyword))            
    titlkeyword();
  else if (CCP4::ccp4_keymatch("LABO", parser->keyword))            
    labokeyword();
  else if (CCP4::ccp4_keymatch("STAT", parser->keyword))            
    statkeyword();
  else if (CCP4::ccp4_keymatch("ALLI", parser->keyword))            
    allikeyword();
  else if (CCP4::ccp4_keymatch("OUTP", parser->keyword))    
    outpkeyword();
  else if (CCP4::ccp4_keymatch("NOHA", parser->keyword))    
    nohakeyword();
  else if (CCP4::ccp4_keymatch("NOUT", parser->keyword))    
    noutkeyword();
  else if (CCP4::ccp4_keymatch("NOUP", parser->keyword))    
    noupkeyword();
  else if (CCP4::ccp4_keymatch("REFA", parser->keyword))  // parameters to refine
    refakeyword();
  else if (CCP4::ccp4_keymatch("MOCC", parser->keyword))  // maximum parameter shifts  
    mocckeyword();                                           
  else if (CCP4::ccp4_keymatch("MXYZ", parser->keyword))            
    mxyzkeyword();
  else if (CCP4::ccp4_keymatch("MBFA", parser->keyword))            
    mbfakeyword();
  else if (CCP4::ccp4_keymatch("MLUZ", parser->keyword))            
    mluzkeyword();
  else if (CCP4::ccp4_keymatch("BETA", parser->keyword))            
    betakeyword();
  else
    return false;
 
 return true;
}

void Input::targkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Input::targkeyword", "Extra text in TARGet keyword will be ignored");
  
  if (!parser->token[1].isstring)	
    Bp3Error("Input::targkeyword", "A string is expected for the TARGet");

  string target = parser->token[1].word; 
  like.Settarget(target);
}

void Input::fillkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Input::fillkeyword", "Extra text in FILL keyword will be ignored");
 
  if (!parser->token[1].isstring)	
    like.Setfill(true);
  else
    if (CCP4::ccp4_keymatch("NO", parser->token[1].word))
      like.Setfill(false);
    else
      like.Setfill(true);
}
  void Input::nodekeyword()
{
  for (int i              = 1; i < parser->ntokens; i++)
  {
    if (CCP4::ccp4_keymatch("CENT", parser->token[i].word))
      if ((i+1)           < parser->ntokens)
      {
	if (parser->token[i+1].isnumber && (parser->token[i+1].value > 0))
	  like.Setcweightsize((unsigned)parser->token[i+1].value);
	else
	  Bp3Error("Input::nodekeyword", "A positive integer is expected for CENT");
      }
      else
	Bp3Error("Input::nodekeyword", "Parsing error");

    if (CCP4::ccp4_keymatch("PHAS", parser->token[i].word))
      if ((i+1)           < parser->ntokens)
      {
	if (parser->token[i+1].isnumber && (parser->token[i+1].value > 0))
	  like.Setpweightsize((unsigned)parser->token[i+1].value);
	else
	  Bp3Error("Input::nodekeyword", "A positive integer is expected for PHAS");
      }
      else
	Bp3Error("Input::nodekeyword", "Parsing error"); 

    if (CCP4::ccp4_keymatch("AMPL", parser->token[i].word))
      if ((i+1)           < parser->ntokens)
      {
	if (parser->token[i+1].isnumber && (parser->token[i+1].value > 0))
	  like.Setafweightsize((unsigned)parser->token[i+1].value);
	else
	  Bp3Error("Input::nodekeyword", "A positive integer is expected for AMPL");
      }
      else
	Bp3Error("Input::nodekeyword", "Parsing error");

    if (CCP4::ccp4_keymatch("SAD", parser->token[i].word))
      if ((i+1)           < parser->ntokens)
      {
	if (parser->token[i+1].isnumber && (parser->token[i+1].value > 0))
	  like.Setsadweightsize((unsigned)parser->token[i+1].value);
	else
	  Bp3Error("Input::nodekeyword", "A positive integer is expected for SAD");
      }
      else
	Bp3Error("Input::nodekeyword", "Parsing error");
  }
}

void Input::outikeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens > 1)
    Bp3Warning("Input::outikeyword", "Extra text in OUTIntensity keyword will be ignored");
  
  like.Setmode("INTENSITY");
}

void Input::outhkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens > 1)
    Bp3Warning("Input::outhkeyword", "Extra text in OUTHkl keyword will be ignored");
  
  like.Setmode("OUTPUTHKL");
}

void Input::outskeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens > 1)
    Bp3Warning("Input::outskeyword", "Extra text in OUTScalepack keyword will be ignored");
  
  like.Setmode("OUTPUTSCALEPACK");
}

void Input::intekeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 1)
    Bp3Warning("Input::intekeyword", "Extra text in INTE keyword will be ignored");

  like.Setinterpolate(true);
}

void Input::verbkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens      == 1)
  {
    xtal.Setverbose(1);
    like.Setverbose(1);
  }
  else if (parser->ntokens == 2)
    if (parser->token[1].isnumber)
    {
      xtal.Setverbose((unsigned)parser->token[1].value);
      like.Setverbose((unsigned)parser->token[1].value);
    }
    else
      Bp3Error("Input::verbkeyword", "A numerical value is expected for VERBose");
  else
    Bp3Warning("Input::verbkeyword", "Extra text in VERBose keyword will be ignored");
}

void Input::warnkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens      == 1)
  {
    xtal.Setwarn(true);
    like.Setwarn(true);
  }
  else
    Bp3Warning("Input::warnkeyword", "Extra text in WARN keyword will be ignored");
}

void Input::titlkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::titlkeyword", "Extra text in TITLe keyword will be ignored");

  if (parser->token[1].isstring)
    like.Settitle(parser->token[1].fullstring);
  else
    Bp3Error("Input::titlkeyword", "A string is expected after the TITLe keyword");
}

void Input::labokeyword()
{

  // first check only checkword
  for (int i                 = 1; i < parser->ntokens; i++)
    if (CCP4::ccp4_keymatch("ONLY",parser->token[i].fullstring))
    {
      like.Setonly(true);
      like.Setallin(false);
    }
  
  for (int i                 = 1; i < parser->ntokens; i++)
  {
    if (CCP4::ccp4_keymatch("F",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzf((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("SIGF",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzsigf((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("FA",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzfa((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("SGFA",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzsigfa((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("EA",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzea((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("SGEA",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzsigea((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("ALPH",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzalpha((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("FB",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzfb((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("PHIB",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzpb((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("HLA",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzhla((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("HLB",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzhlb((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("HLC",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzhlc((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("HLD",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzhld((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("FOM",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzfom((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("FDIFF",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzfdiff((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("PDIFF",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzpdiff((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("FWT",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzfcomb((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if (CCP4::ccp4_keymatch("PHWT",parser->token[i].fullstring))
    {
      if ((i+1)              < parser->ntokens)
	like.Setmtzpcomb((string)(parser->token[i+1].fullstring));
      else
	Bp3Error("Input::labokeyword","Parsing error");
      i++;
    }
    else if(CCP4::ccp4_keymatch("ONLY",parser->token[i].fullstring))
      continue;
    else
      Bp3Error("Input::labokeyword","subkeyword of LABOut not known");
  }
}

void Input::statkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::statkeyword", "Extra text in STAT keyword will be ignored");

  if (parser->token[1].isstring)
    if (CCP4::ccp4_keymatch("OFF", parser->token[1].word))
      like.Setstats(false);
    else if (CCP4::ccp4_keymatch("ON", parser->token[1].word))
      like.Setstats(true);
    else
      Bp3Warning("Input::statkeyword", "The ON or OFF string is expected after STAT");
}

void Input::allikeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 1)
    Bp3Warning("Input::allikeyword", "Extra text in ALLIn keyword will be ignored");
  like.Setallin(true);
}

void Input::outpkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Input::outpkeyword", "Extra text in OUTPut keyword will be ignored");

  if (parser->ntokens  < 2)
    Bp3Error("Input::outpkeyword", "A string is expected after OUTPut keyword");
  else
    like.Setoutput(parser->token[1].fullstring);
}

void Input::nohakeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 1)
    Bp3Warning("Input::nohakeyword", "Extra text in NOHAnd keyword will be ignored");

    like.Sethand(0);
}

void Input::noutkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 1)
    Bp3Warning("Input::noutkeyword", "Extra text in NOUTput keyword will be ignored");

  like.Setnooutput(true);
}

void Input::noupkeyword()
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 1)
    Bp3Warning("Input::noupkeyword", "Extra text in NOUTput keyword will be ignored");

  like.Setupdatesigmah(false);
}

void Input::refakeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 1)
    Bp3Warning("Input::refakeyword", "Extra text in REFAll keyword will be ignored");

  like.Setrefineall(true);
}

void Input::mocckeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Input::mocckeyword", "Extra text in MOCC keyword will be ignored");

  if (!parser->token[1].isnumber && (parser->token[1].value <= ZERO) )
    Bp3Error("Input::mocckeyword", "A positive numerical value is expected for MOCC");

  like.Setmaxocc((double)(parser->token[1].value));
}

void Input::mluzkeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Input::mluzkeyword", "Extra text in MLUZ keyword will be ignored");

  if (!parser->token[1].isnumber && (parser->token[1].value <= ZERO) )
    Bp3Error("Input::mluzkeyword", "A positive numerical value is expected for MLUZ");

  like.Setmaxluzz((double)(parser->token[1].value));
}

void Input::betakeyword() 
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Multicombinput::betakeyword", "Extra text in BETA keyword will be ignored");

  if (parser->ntokens == 2)
    if (parser->token[1].isnumber)
      like.Setbeta((double)parser->token[1].value);
    else
      Bp3Error("Multicombinput::betakeyword", "A numerical value is expected after BETA");
  else
    Bp3Error("Multicombinput::betakeyword", "A numerical value is expected after BETA");

}

void Input::mbfakeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Input::mbfakeyword", "Extra text in MBFA keyword will be ignored");

  if (!parser->token[1].isnumber && (parser->token[1].value <= ZERO) )
    Bp3Error("Input::mbfakeyword", "A positive numerical value is expected for MBFA");

  like.Setmaxbfac((double)(parser->token[1].value));
}

void Input::mxyzkeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Input::mxyzkeyword", "Extra text in MXYZ keyword will be ignored");

  if (!parser->token[1].isnumber && (parser->token[1].value <= ZERO) )
    Bp3Error("Input::mxyzkeyword", "A positive numerical value is expected for MXYZ");

  like.Setmaxxyz((double)(parser->token[1].value));
}
