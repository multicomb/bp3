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
#include <string.h>
#include <time.h>
#include "cmtzlib.h"
#include "crystal.h"
#include "csymlib.h"

bool filter_atom(const clipper::Atom & at)
{
	if ((at.occupancy() == 0) ||
	    (at.u_iso() == 0)
	   ) return true;
	else return false;
}

void  prepare_monomer (clipper::MMonomer &mono_in, clipper::MMonomer &mono_out)
{
	if (mono_in.type() != "HOH") 
		for ( int a = 0; a < mono_in.size(); a++ ) {
    		//std::cout  << mono_in.id() << "\t" << mono_in.type()+"\t"+mono_in[a].id() <<
			//	"\t" << mono_in[a].coord_orth().x() << "\t" << mono_in[a].coord_orth().y() << "\t" << mono_in[a].coord_orth().z() <<
			//	"\t" << mono_in[a].occupancy() << "\t" << mono_in[a].u_iso()  << std::endl;
			if (!(filter_atom(mono_in[a]))) mono_out.insert(mono_in[a]);
		}
}

void  prepare_chain (clipper::MPolymer &mpoly_in, clipper::MPolymer &mpoly_out)
{
   for ( int m = 0; m < mpoly_in.size(); m++ ) {
	clipper::MMonomer mono_tmp;
	mono_tmp.set_id(mpoly_in[m].id());
	mono_tmp.set_type(mpoly_in[m].type());
	mono_tmp.set_seqnum(mpoly_in[m].seqnum());

	prepare_monomer(mpoly_in[m], mono_tmp);

	if (mono_tmp.size() > 0) mpoly_out.insert(mono_tmp);
   }
}

void  prepare_pdb_model (clipper::MiniMol &model_base, clipper::MiniMol &model_base_tmp)
{
  for ( int p = 0; p < model_base_tmp.size(); p++ ) {
	clipper::MPolymer mpoly_tmp;
	mpoly_tmp.set_id(model_base_tmp[p].id());

	prepare_chain(model_base_tmp[p], mpoly_tmp);

	if(mpoly_tmp.size() > 0) model_base.insert(mpoly_tmp);
    //std::cout << model_base[p].id()+"\t"+model_base[p][m].type()+"\t"+model_base[p][m][a].id() << std::endl;
   }
}

// Cell methods

Cell &Cell::Setcell(const float *cellin)
{
  // set cell parameters
  Seta(cellin[0]).Setb(cellin[1]).Setc(cellin[2]);
  
  Setalpha(cellin[3]).Setbeta(cellin[4]).Setgamma(cellin[5]);
  
  Setor2frac();

  Setabcstar();

  return *this;
}

void Cell::Setabcstar()
{  
  work cosalpha = cos(alpha*DEGREEtoRAD);
  work sinalpha = sin(alpha*DEGREEtoRAD);
  work cosbeta  = cos(beta*DEGREEtoRAD);
  work sinbeta  = sin(beta*DEGREEtoRAD);
  work cosgamma = cos(gamma*DEGREEtoRAD);
  work singamma = sin(gamma*DEGREEtoRAD);
  volume        = a*b*c*sqrt(ONE + TWO*cosalpha*cosbeta*cosgamma -
			     cosalpha*cosalpha - cosbeta*cosbeta -
			     cosgamma*cosgamma);
  astar         = b*c*sinalpha/volume;
  bstar         = a*c*sinbeta/volume;
  cstar         = a*b*singamma/volume;
}

void Cell::Setor2frac()
{
  // Sets orthogonal <-> fractional conversion matrix
  or2frac.resize(3);
  frac2or.resize(3);
  for (unsigned i      = 0; i < 3; i++)
  {
    or2frac[i].resize(3, ZERO);
    frac2or[i].resize(3, ZERO);
  }

  work cosalpha        = cos(alpha*DEGREEtoRAD);
  work cosbeta         = cos(beta*DEGREEtoRAD);
  work sinbeta         = sin(beta*DEGREEtoRAD);
  work cosgamma        = cos(gamma*DEGREEtoRAD);
  work singamma        = sin(gamma*DEGREEtoRAD);
  work arg1            = (cosbeta*cosgamma - cosalpha)/(sinbeta*singamma);
  work arg2            = sqrt(ONE - arg1*arg1);

  or2frac[0][0]        = ONE/a;

  if (fabs(cosgamma)   > EPSILON)
    or2frac[0][1]      = -cosgamma/(singamma*a);

  if ( (fabs(cosgamma) > EPSILON) || (fabs(cosbeta) > EPSILON) )
    or2frac[0][2]      = -((cosgamma*sinbeta*arg1+cosbeta*singamma)/
			   (sinbeta*arg2*singamma*a));
  or2frac[1][1]        = ONE/(singamma*b);
  if (fabs(arg1)       > EPSILON)
    or2frac[1][2]      = arg1/(arg2*singamma*b);

  or2frac[2][2]        = ONE/(sinbeta*arg2*c);

  frac2or[0][0]        = a;
  frac2or[0][1]        = cosgamma*b;
  frac2or[0][2]        = cosbeta*c;

  frac2or[1][1]        = singamma*b;
  frac2or[1][2]        = -sinbeta*arg1*c;

  frac2or[2][2]        = sinbeta*arg2*c;  
}

void Cell::Setstolsq(const vector<vector<int> > &miller)
{
  // Calculates (sin(theta)/lambda)^2
  
  for (unsigned r = 0; r < miller.size(); r++)
    stolsq.push_back((pow(or2frac[0][0]*miller[r][0], 2) +
		      pow(or2frac[0][1]*miller[r][0] + 
			  or2frac[1][1]*miller[r][1], 2) +
		      pow(or2frac[0][2]*miller[r][0] + 
			  or2frac[1][2]*miller[r][1] + 
			  or2frac[2][2]*miller[r][2], 2))/FOUR);
}

void Cell::Setstolsq(const vector<vector<int> > &miller, const unsigned r)
{
  // Calculates (sin(theta)/lambda)^2
  stolsq.push_back((pow(or2frac[0][0]*miller[r][0], 2) +
		    pow(or2frac[0][1]*miller[r][0] + 
			or2frac[1][1]*miller[r][1], 2) +
		    pow(or2frac[0][2]*miller[r][0] + 
			or2frac[1][2]*miller[r][1] + 
			or2frac[2][2]*miller[r][2], 2))/FOUR);
}

void Cell::print(const bool verbose) const
{
  printf("Cell Dimensions: %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n\n",a, b, c, alpha, beta, gamma);

  if (verbose)
  {
    printf("Cell volume %8.2f\n",volume);
    printf("Reciprocal space cell edges\n");
    printf("a* = %8.5f  b* = %8.5f  c* = %8.5f \n",astar, bstar, cstar);
  }
  fflush(stdout);
}

// Spacegroup methods

Spacegroup::Spacegroup(const int numb, const char* name) 
{
  CSym::CCP4SPG *spg;
  
  if (numb              > 0)
  {
    number              = numb;
    // do minor checking on parameters
    if (number          > 215)
      Bp3Error("Spacegroup::Spacegroup", "invalid CCP4 Spacegroup number from mtz");
    spg                 = CSym::ccp4spg_load_by_ccp4_num(number);

  }
  else if (name        != NULL)
  {
    spg                 = CSym::ccp4spg_load_by_spgname(name);
    number              = spg->spg_ccp4_num;
  }
  else
    Bp3Error("Spacegroup::Spacegroup", "Space group name or number not given in constructor");
    
  strncpy(&lattice, spg->symbol_old,1);
  checklattice();

  Setsystem();
  int nsym              = spg->nsymop;
  int nsymp             = spg->nsymop_prim;
  hmname                = spg->symbol_xHM;
  pointgroup            = spg->point_group;
  if ((nsym             > 0) && (nsym < 193))
    NSYM                = (unsigned) nsym;
  else
    Bp3Error("Spacegroup::Spacegroup", "invalid NSYM from mtz");

  if ((nsymp            > 0) && (nsymp < 193))
    NSYMP               = (unsigned) nsymp;
  else
    Bp3Error("Spacegroup::Spacegroup", "invalid NSYMP from mtz");

  // allocate and store symmetry matrices
  symrot.resize(NSYM);
  symtran.resize(NSYM);
  
  for (unsigned s       = 0; s < NSYM; s++)
  {
    symrot[s].resize(3);
    symtran[s].resize(3);
    for (int i          = 0; i < 3; i++)
      symrot[s][i].resize(3);  
    for (int i          = 0; i < 3; i++)
    {  
      // storing transpose of rotation matrix
      for (int j        = 0; j < 3; j++)
	symrot[s][i][j] = (int) spg->symop[s].rot[j][i];
      symtran[s][i]     = (double)spg->symop[s].trn[i];
    }
  }
  Setpolar();
  CSym::ccp4spg_free(&spg);
}

Spacegroup &Spacegroup::Setpolar()
{
  // modified from CCP4's symlib.f (that was taken from Alexei Vagin)
  
  for (unsigned i    = 0; i < 3; i++)
    polar[i]         = true;

  double x           = 0.13;
  double y           = 0.17;
  double z           = 0.19;
  
  for (unsigned isym = 0; isym < NSYM; isym++)
  {
    double xx        = (symrot[isym][0][0]*x +
			symrot[isym][0][1]*y +
			symrot[isym][0][2]*z);

    double yy        = (symrot[isym][1][0]*x +
			symrot[isym][1][1]*y +
			symrot[isym][1][2]*z);

    double zz        = (symrot[isym][2][0]*x +
			symrot[isym][2][1]*y +
			symrot[isym][2][2]*z);

    if (fabs(x - xx) > 0.01)
      polar[0]       = false;
    if (fabs(y - yy) > 0.01)
      polar[1]       = false;
    if (fabs(z - zz) > 0.01)
      polar[2]       = false;
  }  
  return *this;
}

bool Spacegroup::specialpos(double &x, double &y, double &z, double &occ) const
{
  double mult(ZERO);
  
  for (unsigned isym   = 0; isym < NSYMP; isym++)
  {
    double xx          = (symrot[isym][0][0]*x +
			  symrot[isym][0][1]*y +
			  symrot[isym][0][2]*z);

    double yy          = (symrot[isym][1][0]*x +
			  symrot[isym][1][1]*y +
			  symrot[isym][1][2]*z);

    double zz          = (symrot[isym][2][0]*x +
			  symrot[isym][2][1]*y +
			  symrot[isym][2][2]*z);

    if ( (fabs(x - xx) < 0.01) && (fabs(y - yy) < 0.01) && (fabs(z - zz) < 0.01) )
    {
      x                = HALF*(x + xx);
      y                = HALF*(y + yy);
      z                = HALF*(z + zz);      
      mult            += ONE;
    }
  }

  if (mult             > ONE)
    occ               /= mult;
  
  return (mult         > ONE);
}

unsigned Spacegroup::enantiomorph() const
{
  // return number of enantiomorph if enantiomorphic else 0
  unsigned enant(0);

  if (number      == 76 )          // P41
    enant          = 78;
  else if (number == 78 )          // P43
    enant          = 76;
  else if (number == 91 )          // P4122
    enant          = 95;
  else if (number == 95 )          // P4322
    enant          = 91;
  else if (number == 92 )          // P41212
    enant          = 96;
  else if (number == 96 )          // P43212
    enant          = 92;
  else if (number == 144)          // P31
    enant          = 145;
  else if (number == 145)          // P32
    enant          = 144;
  else if (number == 151)          // P3112
    enant          = 153;
  else if (number == 153)          // P3212
    enant          = 151;
  else if (number == 152)          // P3121
    enant          = 154;
  else if (number == 154)          // P3221
    enant          = 152;
  else if (number == 169)          // P61
    enant          = 170;
  else if (number == 170)          // P65
    enant          = 169;
  else if (number == 171)          // P62
    enant          = 172;
  else if (number == 172)          // P64
    enant          = 171;
  else if (number == 178)          // P6122
    enant          = 179;
  else if (number == 179)          // P6522
    enant          = 178;
  else if (number == 180)          // P6222
    enant          = 181;
  else if (number == 181)          // P6422
    enant          = 180;
  else if (number == 212)          // P4332
    enant          = 213;
  else if (number == 213)          // P4132
    enant          = 212;
  
  return enant;
}

Spacegroup &Spacegroup::Setsystem()
{
  // set system based on inputted ccp4 spacegroup number

  if ((number        >= 1  ) && (number <= 2  ))
    system            = "Triclinic";
  else if ((number   >= 3  ) && (number <= 15 ))
    system            = "Monoclinic";
  else if ((number   >= 16 ) && (number <= 74 ))
    system            = "Orthorhombic";
  else if ((number   >= 75 ) && (number <= 142))
    system            = "Tetragonal";
  else if ((number   >= 143) && (number <= 167))
    system            = "Trigonal";
  else if ((number   >= 168) && (number <= 194))
    system            = "Hexagonal";
  else if ((number   >= 195) && (number <= 230))
    system            = "Cubic";    
  
  return *this;
}

unsigned Spacegroup::naniso() const
{
  // number of anisotropic scale factors for
  // the given spacegroup
  unsigned pars(0);
  
  if (system      == "Triclinic")
    pars           = 6;
  else if (system == "Monoclinic")
    pars           = 4;
  else if (system == "Orthorhombic")
    pars           = 3;
  else if (system == "Tetragonal")
    pars           = 2;
  else if (system == "Trigonal")
    pars           = 3;
  else if (system == "Hexagonal")
    pars           = 3;
  else if (system == "Cubic")
    pars           = 0;

  return pars;
}

void Spacegroup::checklattice() const
{
  // make sure that the lattice is one of the
  // allowed characters

  if ( !((lattice == 'P') || (lattice == 'A') || 
	 (lattice == 'B') || (lattice == 'C') || 
	 (lattice == 'F') || (lattice == 'I') ||
	 (lattice == 'R') || (lattice == 'H')))
    Bp3Error("Spacegroup::Checklattice","Lattice type not recognized");
}

void Spacegroup::print(const bool verbose) const
{
  printf("Space group information\n");
  printf("Space group name: %s\n", hmname.c_str());
  printf("Space group number: %u\n",number);
  printf("Point group: %s\n",pointgroup.c_str());
  printf("Number of symmetry operators: %u\n",NSYM);
  printf("Number of primitive symmetry operators: %u\n", NSYMP);
  printf("Lattice type: %c\n", lattice);
  printf("Crystal system: %s\n", system.c_str());
  if (polar[0])
    printf("Polar spacegroup: origin not fixed along a axis\n");
  if (polar[1])
    printf("Polar spacegroup: origin not fixed along b axis\n");
  if (polar[2])
    printf("Polar spacegroup: origin not fixed along c axis\n");
  printf("\n");
  fflush(stdout);
}

Crystal::Crystal()
{
  // Crystal defaults
  verbose   = onlycentrics  = onlyacentrics = heavyref = mad = rscale = clipper = invert = warn = false;
  nd        = -1; userpluz.resize(6,false); defaultpluz.resize(6,ONE); pluz.resize(4); refinep.resize(4,false);
  spg_trial = -1;
}

Crystal &Crystal::setwithunknownfiletype(const string type, const vector<string> name)
{
  // name the mtz file MTZIN
  vector<CMtz::MTZ *> MTZIN(sf.size());
  bool allmtz(true), allsca(true);

  printf("Checking if the input reflection file is mtz\n\n");
  for (unsigned d            = 0; d < sf.size(); d++)
  {
    MTZIN[d]                 = CMtz::MtzGet(name[d].c_str(),1);
    if (MTZIN[d])
    {
      sf[d].Setmtzin(name[d].c_str());
      allsca                 = false;
    }
    else
    {
      sf[d].Setscain(name[d].c_str());
      allmtz                 = false;
    }
  }

  if ((!allmtz)             && (!allsca))
    Bp3Error("Crystal::setwithunknownfiletype","both scalepack and mtz files detected");
  else if (allmtz)
    setwithmanymtz(type,false);
  else if (allsca)
  {
    Bp3Warning("Crystal::setwithunknownfiletype","Input reflection file is not mtz, assuming it is a scalepack file.");    
    setwithmanysca(type);
  }
  
  return *this;
}

Crystal &Crystal::setwithmanymtz(const string type, const bool amp)
{
  // name the mtz file MTZIN
  vector<CMtz::MTZ *> MTZIN(sf.size());

  Bp3Warning("Crystal::setwithmanymtz","Assuming the mtz file(s) are all in the CCP4 asymmetric unit");

  for (unsigned d            = 0; d < sf.size(); d++)
  {
    MTZIN[d]                 = CMtz::MtzGet(sf[d].mtzin.c_str(),1);
    if (!MTZIN[d])
    {
      const string temp("MTZ " + sf[d].mtzin + " not found");
      Bp3Error("Crystal::setwithmanymtz",temp);
    }
  }

  // check if spacegroups match
  for (unsigned d1           = 0   ; d1 < sf.size(); d1++)
    for (unsigned d2         = d1+1; d2 < sf.size(); d2++)
      if (CMtz::MtzSpacegroupNumber(MTZIN[d1]) != CMtz::MtzSpacegroupNumber(MTZIN[d2]))
      {
	const string temp("Spacegroups do not match for mtz files " + sf[d1].mtzin + " and " + sf[d2].mtzin);
	Bp3Error("Crystal::setwithmanymtz",temp);
      }
  
  // Set up spacegroup
  CMtz::SYMGRP mtzsym        = MTZIN[0]->mtzsymm;
  sg                         = Spacegroup(0,mtzsym.spcgrpname);

  vector<unsigned> datapcol(sf.size(),0);
  vector<unsigned> devpcol(sf.size(),0);
  vector<unsigned> datamcol(sf.size(),0);
  vector<unsigned> devmcol(sf.size(),0);
  vector<vector<unsigned> > modelfcol(sf.size());
  vector<vector<unsigned> > modelpcol(sf.size());

  vector<vector<CMtz::MTZCOL * > > col(sf.size());
  vector<CMtz::MTZSET * > set(sf.size());
  vector<CMtz::MTZXTAL * > xtl(sf.size());
  vector< vector<int> > xtalnumber(sf.size());

  for (unsigned d            = 0; d < sf.size(); d++)
  {
    modelfcol[d].resize(sf[d].sfmodel.size(),0);    
    modelpcol[d].resize(sf[d].sfmodel.size(),0);
    for (int i               = 0; i < CMtz::MtzNxtal(MTZIN[d]); i++) 
    {
      xtl[d]                 = CMtz::MtzIxtal(MTZIN[d],i);
      for (int j             = 0; j < CMtz::MtzNsetsInXtal(xtl[d]);  j++) 
      {
	set[d]               = CMtz::MtzIsetInXtal(xtl[d],j);
	for (int k           = 0; k < CMtz::MtzNcolsInSet(set[d]); k++)
	{
	  col[d].push_back(CMtz::MtzIcolInSet(set[d],k));
	  xtalnumber[d].push_back(i);
	}
      }
    }
  }


  // prefer intensities over ampltitudes - check if they exist.
 
  bool ampl(amp);

  if (!ampl)
  {
    for (unsigned d          = 0; d < sf.size(); d++)
      for (unsigned j        = 0; j < col[d].size(); j++)
	if (!sf[d].sdatap.size())
	  if ( (strstr(col[d][j]->type,"G") != NULL) && (strstr(col[d][j]->label,"+") != NULL) )
	    if ( (sf[d].type == "AMPLITUDE") || (!sf[d].type.size()) )
	    {
	      ampl           = true;
	      break;
	    }
    if (!ampl) 
      for (unsigned d          = 0; d < sf.size(); d++)
        for (unsigned j        = 0; j < col[d].size(); j++)
	  if (!sf[d].sdatap.size())
	    if ( (strstr(col[d][j]->type,"K") != NULL) && (strstr(col[d][j]->label,"+") != NULL) && !amp)
	    {
	      for (unsigned d1 = 0; d1 < sf.size(); d1++)
	        sf[d1].type    = "INTENSITY";
	      break;
	    }
  }
  
  // see if the inputted column set are present and set there column numbers
  // if no column labels are given, try to set it if only one type is given

  vector<unsigned> iplus(sf.size(),0), iminus(sf.size(),0), fplus(sf.size(),0), fminus(sf.size(),0);
  vector<unsigned> siplus(sf.size(),0), siminus(sf.size(),0), sfplus(sf.size(),0), sfminus(sf.size(),0);
  
  for (unsigned d            = 0; d < sf.size(); d++)
    for (unsigned j          = 0; j < col[d].size(); j++)
    {
      if (sf[d].sdatap.size())
      {
	if (!strcmp(col[d][j]->label, sf[d].sdatap.c_str()))
	  datapcol[d]        = j;
      }
      else if ( (strstr(col[d][j]->type,"K") != NULL) && (strstr(col[d][j]->label,"+") != NULL) && !ampl)
      {
	iplus[d]++;
	datapcol[d]          = j;
	sf[d].sdatap         = col[d][j]->label;
	printf("Setting I+ for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdatap.c_str());
	for (unsigned d      = 0; d < sf.size(); d++)
	  sf[d].type         = "INTENSITY";

      }
      else if ( (strstr(col[d][j]->type,"G") != NULL) && (strstr(col[d][j]->label,"+") != NULL) && !iplus[d] && (sf[d].type != "INTENSITY"))
      {
	fplus[d]++;
	datapcol[d]          = j;
	sf[d].sdatap         = col[d][j]->label;
	printf("Setting F+ for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdatap.c_str());
	for (unsigned d      = 0; d < sf.size(); d++)
	  sf[d].type         = "AMPLITUDE";
      }
      
      if (sf[d].sdevp.size())
      {
	if (!strcmp(col[d][j]->label, sf[d].sdevp.c_str()))
	  devpcol[d]         = j;
      }
      else if ( (strstr(col[d][j]->type,"M") != NULL) && (strstr(col[d][j]->label,"+") != NULL) && !ampl)
      {
	siplus[d]++;
	devpcol[d]           = j;
	sf[d].sdevp          = col[d][j]->label;
	printf("Setting SIGI+ for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdevp.c_str());

      }
      else if ( (strstr(col[d][j]->type,"L") != NULL) && (strstr(col[d][j]->label,"+") != NULL) && !siplus[d] && (sf[d].type != "INTENSITY"))
      {
	sfplus[d]++;
	devpcol[d]           = j;
	sf[d].sdevp          = col[d][j]->label;
	printf("Setting SIGF+ for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdevp.c_str());
      }

      if (sf[d].sdatam.size())
      {
	if (!strcmp(col[d][j]->label, sf[d].sdatam.c_str()))
	  datamcol[d]        = j;
      }
      else if ( (strstr(col[d][j]->type, "K") != NULL) && (strstr(col[d][j]->label,"-") != NULL) && !ampl)
      {
	iminus[d]++;
	datamcol[d]          = j;
	sf[d].sdatam         = col[d][j]->label;
	printf("Setting I- for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdatam.c_str());
	sf[d].anomalous      = true;
      }
      else if ( (strstr(col[d][j]->type, "G") != NULL) && (strstr(col[d][j]->label,"-") != NULL) && !iminus[d] && (sf[d].type != "INTENSITY"))
      {
	fminus[d]++;
	datamcol[d]          = j;
	sf[d].sdatam         = col[d][j]->label;
	printf("Setting F- for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdatam.c_str());
	sf[d].anomalous      = true;
      }

      if (sf[d].sdevm.size())
      {
	if (!strcmp(col[d][j]->label, sf[d].sdevm.c_str()))
	  devmcol[d]         = j;
      }
      else if ( (strstr(col[d][j]->type, "M") != NULL) && (strstr(col[d][j]->label,"-") != NULL) && iminus[d])
      {
	siminus[d]++;
	devmcol[d]           = j;
	sf[d].sdevm          = col[d][j]->label;
	printf("Setting SIGI- for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdevm.c_str());
      }
      else if ( (strstr(col[d][j]->type, "L") != NULL) && (strstr(col[d][j]->label,"-") != NULL) && !siminus[d] && (sf[d].type != "INTENSITY"))
      {
	sfminus[d]++;
	devmcol[d]           = j;
	sf[d].sdevm          = col[d][j]->label;
	printf("Setting SIGF- for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdevm.c_str());
      }

      for (unsigned m        = 0; m < sf[d].sfmodel.size(); m++)
      {
	if (sf[d].sfmodel[m].size())
	  if (!strcmp(col[d][j]->label, sf[d].sfmodel[m].c_str()))
	    modelfcol[d][m]  = j;

	if (sf[d].spmodel[m].size())
	  if (!strcmp(col[d][j]->label, sf[d].spmodel[m].c_str()))
	    modelpcol[d][m]   = j;
      }
    }

  for (unsigned d            = 0; d < sf.size(); d++)
  {
    // check if there is a non-anomalous data set
    if (!datapcol[d]        && !sf[d].sdatap.size())
      for (unsigned j        = 0; j < col[d].size(); j++)
	if (!amp && (strstr(col[d][j]->type,"J") != NULL))
	{
	  iplus[d]++;
	  datapcol[d]        = j;
	  sf[d].sdatap       = col[d][j]->label;
	  sf[d].type         = "INTENSITY";
	  printf("Setting I for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdatap.c_str());
	}
	else if ((strstr(col[d][j]->type,"F") != NULL) && !iplus[d])
	{
	  fplus[d]++;
	  datapcol[d]        = j;
	  sf[d].sdatap       = col[d][j]->label;
	  printf("Setting F for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdatap.c_str());
	  sf[d].type         = "AMPLITUDE";
	}
	else if ((strstr(col[d][j]->type,"Q") != NULL) && iplus[d])
	{
	  siplus[d]++;
	  devpcol[d]         = j;
	  sf[d].sdevp        = col[d][j]->label;
	  printf("Setting SIGI for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdevp.c_str());
	}
	else if ((strstr(col[d][j]->type,"Q") != NULL) && fplus[d])
	{
	  sfplus[d]++;
	  devpcol[d]         = j;
	  sf[d].sdevp        = col[d][j]->label;
	  printf("Setting SIGF for data set %s to %s\n\n",sf[d].name.c_str(),sf[d].sdevp.c_str());
	}
    
    if ( (iplus[d]           > 1) || (iminus[d] > 1) || (siplus[d] > 1) || (siminus[d] > 1) )
      Bp3Error("Crystal::setwithmanymtz", "More than one I+, I-, SIGI+ or SIGI- minus found - please use COLUmn keyword");

    if ( (fplus[d]           > 1) || (fminus[d] > 1) || (sfplus[d] > 1) || (sfminus[d] > 1) )
      Bp3Error("Crystal::setwithmanymtz", "More than one F+, F-, SIGF+ or SIGF- minus found - please use COLUmn keyword");

    if ( (fplus[d]           > 1) || (fminus[d] > 1) || (sfplus[d] > 1) || (sfminus[d] > 1) )
      Bp3Error("Crystal::setwithmanymtz", "More than one F+, F-, SIGF+ or SIGF- minus found - please use COLUmn keyword");
  }
  
  // now, check to see that the mtz file contains the label that the user asked for and only one type is given
  string t(sf[0].type);

  for (unsigned d            = 0; d < sf.size(); d++)
  {    
    if (!datapcol[d])
      Bp3Error("Crystal::setwithmanymtz", "Column for F/F+ or I/I+:" + sf[d].sdatap + " was not found in mtz");
    if (!devpcol[d])
      Bp3Error("Crystal::setwithmanymtz", "Column for SF/SF+ or SI/SI+:" + sf[d].sdevp + " was not found in mtz");
    if (sf[d].anomalous)
    {
      if (!datamcol[d])
	Bp3Error("Crystal::setwithmanymtz", "Column for F- or I-:" + sf[d].sdatam + " was not found in mtz");
      if (!devmcol[d])
	Bp3Error("Crystal::setwithmanymtz", "Column for SF- or SI-:" + sf[d].sdevm + " was not found in mtz");
    }
    if (sf[d].type         != t)
      Bp3Error("Crystal::setwithmanymtz", "Give only amplitudes or intensities, not both");
    
    if (invert)
    {
      string tempda         = sf[d].Getsdatap();
      string tempde         = sf[d].Getsdevp();
      sf[d].Setsdatap(sf[d].Getsdatam());
      sf[d].Setsdevp(sf[d].Getsdevm());
      sf[d].Setsdatam(tempda);
      sf[d].Setsdevm(tempde);
    }  
  }
  
  // maximum number of reflections in mtz file
  unsigned md(0);
  maxref                     = CMtz::MtzNref(MTZIN[0]);

  vector<int> ind_xtal(sf.size()), ind_set(sf.size());
  vector<vector<int> > ind_col(sf.size());
  
  for (unsigned d            = 0; d < sf.size(); d++)
  {
    ind_col[d].resize(3);
    if (maxref               < (unsigned)(CMtz::MtzNref(MTZIN[d])))
    {
      maxref                 = (unsigned)(CMtz::MtzNref(MTZIN[d]));
      md                     = d;
    }
    if (!CMtz::MtzFindInd(MTZIN[d], &ind_xtal[d], &ind_set[d], &ind_col[d][0]))
      Bp3Error("Crystal::setwithmanymtz", "HKL could not be found in mtz");
  }
  
  maxselref                  = maxref;
  resize();

  Setdataset();

  // setup the different cells
  for (unsigned c            = 0; c < cell.size(); c++)
    if (cell[c].a           == 0.0)
    { 	
      unsigned d(dataset[c][0]);
      for (unsigned j        = 0; j < col[d].size(); j++)
	if (!strcmp(col[d][j]->label, sf[d].sdatap.c_str()))
        {
	  xtl[d]             = CMtz::MtzIxtal(MTZIN[d],xtalnumber[d][j]);
	  cell[c].Setcell(xtl[d]->cell);
	}  	 
    }

  vector<work> highestres(sf.size(), 1000.0);
  vector<work> lowestres(sf.size(), ZERO);
  
  for (unsigned d            = 0; d < sf.size(); d++)
    sf[d].nref               = 0;

  for (unsigned r            = 0; r < maxref; r++)
  {
    miller[r][0]             = (int)MTZIN[md]->xtal[ind_xtal[md]]->set[ind_set[md]]->col[ind_col[md][0]]->ref[r];
    miller[r][1]             = (int)MTZIN[md]->xtal[ind_xtal[md]]->set[ind_set[md]]->col[ind_col[md][1]]->ref[r];
    miller[r][2]             = (int)MTZIN[md]->xtal[ind_xtal[md]]->set[ind_set[md]]->col[ind_col[md][2]]->ref[r];

    // set cell
    for (unsigned c          = 0; c < cell.size(); c++)
      cell[c].Setstolsq(miller,r);

    // set epsilon and centricity
    Setepsilon(r).Setcentric(r);

    if (Sysabs(r))
    {
      if (verbose)
	Bp3Warning("Crystal::setwithmanymtz","Systematically absent reflection found - it won't be used or phased");
      for (unsigned d        = 0; d < sf.size(); d++)
      {
	sf[d].datap[r]       = NOTUSED;
 	sf[d].devp[r]        = NOTUSED;
	if (sf[d].anomalous)
	{
	  sf[d].datam[r]     = NOTUSED;
	  sf[d].devm[r]      = NOTUSED;
	}
      }
    }
    else
      for (unsigned d        = 0; d < sf.size(); d++)
      {
	unsigned r2(0);

	if (d               != md)               
	  while (r2          < (unsigned)(CMtz::MtzNref(MTZIN[d])))
	    if ((miller[r][0] == (int)MTZIN[d]->xtal[ind_xtal[d]]->set[ind_set[d]]->col[ind_col[d][0]]->ref[r2]) &&
		(miller[r][1] == (int)MTZIN[d]->xtal[ind_xtal[d]]->set[ind_set[d]]->col[ind_col[d][1]]->ref[r2]) &&
		(miller[r][2] == (int)MTZIN[d]->xtal[ind_xtal[d]]->set[ind_set[d]]->col[ind_col[d][2]]->ref[r2]))
	      break;
	    else
	      r2++;
	else
	  r2                 = r;
	
	if (r2              >= (unsigned)(CMtz::MtzNref(MTZIN[d]) - 1))
	  break;
	
	if (!CMtz::ccp4_ismnf(MTZIN[d],col[d][datapcol[d]]->ref[r2]) && (col[d][devpcol[d]]->ref[r2] > ZERO) )
	{
	  // do not store reflections with non-positive sigmas
	  sf[d].datap[r]     = (work)col[d][datapcol[d]]->ref[r2];
	  sf[d].devp[r]      = (work)col[d][devpcol[d]]->ref[r2];
	}
	else
	{
	  sf[d].datap[r]     = NOTUSED;
	  sf[d].devp[r]      = NOTUSED;
	}

	for (unsigned m           = 0; m < sf[d].sfmodel.size(); m++)
	  if (sf[d].sfmodel[m].size())
	    if (!CMtz::ccp4_ismnf(MTZIN[d],col[d][modelfcol[d][m]]->ref[r2]) &&
		!CMtz::ccp4_ismnf(MTZIN[d],col[d][modelpcol[d][m]]->ref[r2])   )
	    {
	      sf[d].fmodel[r][m]  = (work)col[d][modelfcol[d][m]]->ref[r2];
	      sf[d].pmodel[r][m]  = (work)col[d][modelpcol[d][m]]->ref[r2];
	    }
	    else
	    {
	      sf[d].fmodel[r][m]  = NOTUSED;
	      sf[d].pmodel[r][m]  = NOTUSED;
	    }
	  
	
	if (sf[d].anomalous)
	  if (!CMtz::ccp4_ismnf(MTZIN[d],col[d][datamcol[d]]->ref[r2]) && (col[d][devmcol[d]]->ref[r2] > ZERO) ) 
	  {
	    sf[d].datam[r]   = (work)col[d][datamcol[d]]->ref[r2];
	    sf[d].devm[r ]   = (work)col[d][devmcol[d]]->ref[r2];
	  }
	  else if (centric[r] && (sf[d].datap[r] > NOTUSED) )  
	  {
	    // ***NSP - default in many programs is to have F- == absent for centrics
	    sf[d].datam[r]   = sf[d].datap[r];
	    sf[d].devm[r]    = sf[d].devp[r];
	  }
	  else
	  {
	    sf[d].datam[r]   = NOTUSED;
	    sf[d].devm[r]    = NOTUSED;
	  }
      }

    // sigma cut off's and resolution limits

    for (unsigned d          = 0; d < sf.size(); d++)
    {
      if (sf[d].datap[r]     > NOTUSED)
      {
	if (sf[d].datap[r]/sf[d].devp[r] < sf[d].sigmacut)
	  sf[d].datap[r]     = sf[d].devp[r] = NOTUSED;
	if ( (d              > 0) && sf[0].use(r) )
	  if (fabs(sf[d].datamean(r) - sf[0].datamean(r))/
	      sqrt(sf[d].devmean(r)*sf[d].devmean(r) +
		   sf[0].devmean(r)*sf[0].devmean(r)) < sf[d].isosigmacut)
	  {
	    sf[d].datap[r]   = sf[d].devp[r] = NOTUSED;
	    if (sf[d].anomalous)
	      sf[d].datam[r] = sf[d].devm[r] = NOTUSED;
	  }
	
      	/*
	if (sf[d].anomalous)
	{
	  if (sf[d].datam[r] > NOTUSED)
	    if (sf[d].datam[r]/sf[d].devm[r] < sf[d].sigmacut)
	      sf[d].datam[r] = sf[d].devm[r] = NOTUSED;
	  if (sf[d].anouse(r) )
	    if (fabs(sf[d].dano(r))/sf[d].sigdano(r) < sf[d].anosigmacut)
	      sf[d].datam[r] = sf[d].devm[r] = sf[d].datap[r] = sf[d].devp[r] = NOTUSED;
	}
	*/
      }
      work resol             = Getres(d,r);
      // calculate the highest and lowest resolution limits
      if (sf[d].use(r))
      {
 	highestres[d]        = std::min(resol, highestres[d]);
	lowestres[d]         = std::max(resol, lowestres[d]);
      }

      // if the user has specified resolution limits,
      // remove any reflections outside of the boundary 
      if (sf[d].hires       != sf[d].MAXRES)
	if (resol            < sf[d].hires)
	{
	  sf[d].datap[r]     = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]   = NOTUSED;
	}
      if (sf[d].lowres      != sf[d].MINRES)
	if (resol            > sf[d].lowres)
	{
	  sf[d].datap[r]     = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]   = NOTUSED;
	}
      if (sf[d].use(r))
	sf[d].nref++;
    }
  }

  for (unsigned d            = 0; d < sf.size(); d++)
  {
    // user hasn't set resolution limits, so set it here.
    if (sf[d].hires         == sf[d].MAXRES)
      sf[d].hires            = highestres[d];
    if (sf[d].lowres        == sf[d].MINRES)
      sf[d].lowres           = lowestres[d];
  
    if (MTZIN[d]) 
      CMtz::MtzFree(MTZIN[d]);  

    if (sf[d].nwave          > 0)
      mad                    = true;
  
    // set nbins, if not inputted by user
    if (!sf[d].nbins)        
      sf[d].Setnbins();
    sf[d].resizebin();
    sf[d].Setshell();
  }

  if (mad)     // set nbins == 5 (if user hasn't set it)
    Setonlyacentrics(true);  // ***NSP

  Setselected(type);

  for (unsigned d            = 0; d < sf.size(); d++)
    sf[d].checkparameters();

  return *this;
}

Crystal &Crystal::setwithmanysca(const string type)
{
  Bp3Warning("Crystal::setwithmanysca","Assuming the sca file(s) are all in the CCP4 asymmetric unit");
  vector<FILE *> scafile(sf.size(),NULL);
  char line[71];
  vector<string> sgrp(sf.size());
  vector<vector<float> > incell(sf.size());
  
  for (unsigned d            = 0; d < sf.size(); d++)
  {
    incell[d].resize(6,ZERO);
    
    scafile[d]               = fopen(sf[d].scain.c_str(),"r");
    if (scafile[d]          == NULL)
    {
      const string temp("Scalepack file " + sf[d].scain + " doesn't exist");
      Bp3Error("Crystal::setwithmanysca",temp);  
    }
    fgets(line,70,scafile[d]);
    fgets(line,70,scafile[d]);
    fgets(line,70,scafile[d]);
    char spacegroup[25];
    sscanf(line," %9f %9f %9f %9f %9f %9f %s",&incell[d][0],&incell[d][1],&incell[d][2],&incell[d][3],&incell[d][4],&incell[d][5],&sgrp[d][0]);
    for (unsigned c          = 0; c < 6; c++)
      if (incell[d][c]      <= ZERO)
	Bp3Error("Crystal::setwithmanysca","The cell dimensions from the scalepack file were not read correctly");
    
    sf[d].type              == "INTENSITY";
    
  }
  
  // check if spacegroups match
  for (unsigned d1           = 0   ; d1 < sf.size(); d1++)
    for (unsigned d2         = d1+1; d2 < sf.size(); d2++)
      if (sgrp[d1]          != sgrp[d2])
      {
	const string temp("Spacegroups do not match for mtz files " + sf[d1].scain + " and " + sf[d2].scain);
	Bp3Error("Crystal::setwithmanysca",temp);
      }
  
  // Set up spacegroup
  sg                         = Spacegroup(-1,sgrp[0].c_str());

  // maximum number of reflections

  unsigned md(0);
  vector<unsigned> nline(sf.size(),ZERO);
  for (unsigned d            = 0; d < sf.size()   ; d++)
  {
    while (!feof(scafile[d]))
    {
      fgets(line,70,scafile[d]);
      if (strlen(line)       > 30)
	sf[d].anomalous      = true;
      nline[d]++;
    }
    fclose(scafile[d]);
    // set labels
    sf[d].type               = "INTENSITY";
    string i("I"), sigi("SIGI");
    char datanumber[3];
    sprintf(datanumber,"%1d",d+1);
    if (sf[d].anomalous)
    {
      sf[d].Setsdatap(i + datanumber + "(+)");
      sf[d].Setsdatam(i + datanumber + "(-)");
      sf[d].Setsdevp(sigi + datanumber + "(+)");
      sf[d].Setsdevm(sigi + datanumber + "(-)");
    }
    else
    {
      sf[d].Setsdatap(i + datanumber);
      sf[d].Setsdevp(sigi + datanumber);
    }
  }
  
  maxref                     = nline[0];
  for (unsigned d            = 0  ; d  < sf.size()  ; d++)
    for (unsigned d1         = d+1; d1 < sf.size()  ; d1++)
      if (nline[d1]          > nline[d])
      {
	maxref               = nline[d1];
	md                   = d1;
      }
  
  maxselref                  = maxref;
  resize();

  Setdataset();
  
  // setup the different cells
  for (unsigned c            = 0; c < cell.size(); c++)
    if (cell[c].a           == 0.0)
    { 	
      unsigned d(dataset[c][0]);
      cell[c].Setcell(&incell[d][0]);
    }
  
  vector<work> highestres(sf.size(), 1000.0);
  vector<work> lowestres(sf.size(), ZERO);
  
  for (unsigned d            = 0; d < sf.size(); d++)
    sf[d].nref               = 0;

  for (unsigned d            = 0; d < sf.size(); d++)
  {
    scafile[d]               = fopen(sf[d].scain.c_str(),"r");
    fgets(line,70,scafile[d]);
    fgets(line,70,scafile[d]);
    fgets(line,70,scafile[d]);
  }
  
  for (unsigned r            = 0; r < maxref; r++)
  {
    vector<int> h(sf.size()),k(sf.size()),l(sf.size());
    vector<work> ip(sf.size(), NOTUSED), sigip(sf.size(),NOTUSED), im(sf.size(),NOTUSED), sigim(sf.size(),NOTUSED);

    for (unsigned d          = 0; d < sf.size(); d++)
      if (fgets(line,70,scafile[d]) != NULL)
	sscanf(line,"%4d%4d%4d%8f%8f%8f%8f",&h[d],&k[d],&l[d],&ip[d],&sigip[d],&im[d],&sigim[d]);

    miller[r][0]             = h[md];
    miller[r][1]             = k[md];
    miller[r][2]             = l[md];

    // set cell
    for (unsigned c            = 0; c < cell.size(); c++)
      cell[c].Setstolsq(miller,r);

    // set epsilon and centricity
    Setepsilon(r).Setcentric(r);

    if (Sysabs(r))
    {
      if (verbose)
	Bp3Warning("Crystal::setwithmanysca","Systematically absent reflection found - it won't be used or phased");
      for (unsigned d        = 0; d < sf.size(); d++)
      {
	sf[d].datap[r]       = NOTUSED;
 	sf[d].devp[r]        = NOTUSED;
	if (sf[d].anomalous)
	{
	  sf[d].datam[r]     = NOTUSED;
	  sf[d].devm[r]      = NOTUSED;
	}
      }
    }
    else
      for (unsigned d        = 0; d < sf.size(); d++)
      {
	unsigned r2(0);

	if (d               != md)               
	  while (r2          < maxref)
	    if ((miller[r][0] == h[d]) &&
		(miller[r][1] == k[d]) &&
		(miller[r][2] == l[d]))
	      break;
	    else
	      r2++;
	else
	  r2                 = r;
	
	if (r2              >= (unsigned)(maxref - 1))
	  break;
	
	if ( (ip[d]          > NOTUSED) && (sigip[d] > ZERO) )
 	{
	  // do not store reflections with non-positive sigmas
	  sf[d].datap[r]     = (work)ip[d];
	  sf[d].devp[r]      = (work)sigip[d];
	}
	else
	{
	  sf[d].datap[r]     = NOTUSED;
	  sf[d].devp[r]      = NOTUSED;
	}
  
	if (sf[d].anomalous)
	  if ( (im[d]        > NOTUSED) && (sigim[d] > ZERO) ) 
	  {
	    sf[d].datam[r]   = (work)im[d];
	    sf[d].devm[r ]   = (work)sigim[d];
	  }
	  else if (centric[r] && (sf[d].datap[r] > NOTUSED) )  
	  {
	    // ***NSP - default in many programs is to have F- == absent for centrics
	    sf[d].datam[r]   = sf[d].datap[r];
	    sf[d].devm[r]    = sf[d].devp[r];
	  }
	  else
	  {
	    sf[d].datam[r]   = NOTUSED;
	    sf[d].devm[r]    = NOTUSED;
	  }
      }

    // sigma cut off's and resolution limits

    for (unsigned d          = 0; d < sf.size(); d++)
    {      
      if (sf[d].datap[r]     > NOTUSED)
      {
	if (sf[d].datap[r]/sf[d].devp[r] < sf[d].sigmacut)
	  sf[d].datap[r]     = sf[d].devp[r] = NOTUSED;
	if ( (d              > 0) && sf[0].use(r) )
	  if (fabs(sf[d].datamean(r) - sf[0].datamean(r))/
	      sqrt(sf[d].devmean(r)*sf[d].devmean(r) +
		   sf[0].devmean(r)*sf[0].devmean(r)) < sf[d].isosigmacut)
	  {
	    sf[d].datap[r]   = sf[d].devp[r] = NOTUSED;
	    if (sf[d].anomalous)
	      sf[d].datam[r] = sf[d].devm[r] = NOTUSED;
	  }
	
      	/*
	if (sf[d].anomalous)
	{
	  if (sf[d].datam[r] > NOTUSED)
	    if (sf[d].datam[r]/sf[d].devm[r] < sf[d].sigmacut)
	      sf[d].datam[r] = sf[d].devm[r] = NOTUSED;
	  if (sf[d].anouse(r) )
	    if (fabs(sf[d].dano(r))/sf[d].sigdano(r) < sf[d].anosigmacut)
	      sf[d].datam[r] = sf[d].devm[r] = sf[d].datap[r] = sf[d].devp[r] = NOTUSED;
	}
	*/
      }
      work resol             = Getres(d,r);
      // calculate the highest and lowest resolution limits
      if (sf[d].use(r))
      {
 	highestres[d]        = std::min(resol, highestres[d]);
	lowestres[d]         = std::max(resol, lowestres[d]);
      }

      // if the user has specified resolution limits,
      // remove any reflections outside of the boundary 
      if (sf[d].hires       != sf[d].MAXRES)
	if (resol            < sf[d].hires)
	{
	  sf[d].datap[r]     = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]   = NOTUSED;
	}
      if (sf[d].lowres      != sf[d].MINRES)
	if (resol            > sf[d].lowres)
	{
	  sf[d].datap[r]     = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]   = NOTUSED;
	}
      if (sf[d].use(r))
	sf[d].nref++;
    }
  }

  for (unsigned d            = 0; d < sf.size(); d++)
  {
    // user hasn't set resolution limits, so set it here.
    if (sf[d].hires         == sf[d].MAXRES)
      sf[d].hires            = highestres[d];
    if (sf[d].lowres        == sf[d].MINRES)
      sf[d].lowres           = lowestres[d];
  
    fclose(scafile[d]);
    
    if (sf[d].nwave          > 0)
      mad                    = true;
  
    // set nbins, if not inputted by user
    if (!sf[d].nbins)        
      sf[d].Setnbins();
    sf[d].resizebin();
    sf[d].Setshell();
  }

  if (mad)     // set nbins == 5 (if user hasn't set it)
    Setonlyacentrics(true);  // ***NSP

  Setselected(type);

  for (unsigned d            = 0; d < sf.size(); d++)
    sf[d].checkparameters();

  return *this;
}

Crystal &Crystal::setwithmtz(const string type)
{
  // name the mtz file MTZIN
  CMtz::MTZ *MTZIN               = CMtz::MtzGet("HKLIN",1);

  if (!MTZIN)
    Bp3Error("Crystal::setwithmtz","HKLIN not found");

  // Set up spacegroup
  CMtz::SYMGRP mtzsym            = MTZIN->mtzsymm;
  sg                             = Spacegroup(0,mtzsym.spcgrpname);  

  vector<unsigned> datapcol(sf.size(),0);
  vector<unsigned> devpcol(sf.size(),0);
  vector<unsigned> datamcol(sf.size(),0);
  vector<unsigned> devmcol(sf.size(),0);
  vector<vector<unsigned> > modelfcol(sf.size());
  vector<vector<unsigned> > modelpcol(sf.size());
  for (unsigned d                = 0; d < sf.size(); d++)
    if (sf[d].sfmodel.size())
    {
      modelfcol[d].resize(sf[d].sfmodel.size(),0);
      modelpcol[d].resize(sf[d].sfmodel.size(),0);    
    }
  vector<unsigned> phicol(sf.size(),0);
  vector<unsigned> fomcol(sf.size(),0);
  vector<unsigned> hlacol(sf.size(),0);
  vector<unsigned> hlbcol(sf.size(),0);
  vector<unsigned> hlccol(sf.size(),0);
  vector<unsigned> hldcol(sf.size(),0);

  vector<CMtz::MTZCOL * > col;
  CMtz::MTZSET *set;
  CMtz::MTZXTAL *xtl;
  vector<int> xtalnumber;

  for (int i                     = 0; i < CMtz::MtzNxtal(MTZIN); i++) 
  {
    xtl                          = CMtz::MtzIxtal(MTZIN,i);
    for (int j                   = 0; j < CMtz::MtzNsetsInXtal(xtl);  j++) 
    {
      set                        = CMtz::MtzIsetInXtal(xtl,j);
      for (int k                 = 0; k < CMtz::MtzNcolsInSet(set); k++)
      {
	col.push_back(CMtz::MtzIcolInSet(set,k));
	xtalnumber.push_back(i);
      }
    }
  }

  // see if the inputted column set are present and set there column numbers
  for (unsigned j                = 0; j < col.size(); j++)
  {
    for (unsigned d              = 0; d < sf.size(); d++)
      if (!strcmp(col[j]->label, sf[d].sdatap.c_str()))
	datapcol[d]              = j;
      else if (!strcmp(col[j]->label, sf[d].sdevp.c_str()))
	devpcol[d]               = j;
      else if (!strcmp(col[j]->label, sf[d].sfom.c_str()))
	fomcol[d]                = j;
      else if (!strcmp(col[j]->label, sf[d].sphi.c_str()))
	phicol[d]                = j;
      else if (!strcmp(col[j]->label, sf[d].shla.c_str()))
	hlacol[d]                = j;
      else if (!strcmp(col[j]->label, sf[d].shlb.c_str()))
	hlbcol[d]                = j;
      else if (!strcmp(col[j]->label, sf[d].shlc.c_str()))
	hlccol[d]                = j;
      else if (!strcmp(col[j]->label, sf[d].shld.c_str()))
	hldcol[d]                = j;
      else if (!strcmp(col[j]->label, sf[d].sdatam.c_str()))
	datamcol[d]              = j;
      else if (!strcmp(col[j]->label, sf[d].sdevm.c_str()))
	devmcol[d]               = j;
      else
	for (unsigned m          = 0; m < sf[d].sfmodel.size(); m++)
	  if (!strcmp(col[j]->label, sf[d].sfmodel[m].c_str()))
	    modelfcol[d][m]      = j;
 	  else if (!strcmp(col[j]->label, sf[d].spmodel[m].c_str()))
	    modelpcol[d][m]      = j;
  }

  // now, check to see that the mtz file contains the label that the user asked for

  for (unsigned d                = 0; d < sf.size(); d++)
  {    
    if (!datapcol[d])
      Bp3Error("Crystal::setwithmtz", "Column for F/F+ or I/I+:" + sf[d].sdatap + " was not found in mtz");
    if (!devpcol[d])
      Bp3Error("Crystal::setwithmtz", "Column for SF/SF+ or SI/SI+:" + sf[d].sdevp + " was not found in mtz");
    if (sf[d].datam.size())
    {
      if (!datamcol[d])
	Bp3Error("Crystal::setwithmtz", "Column for F- or I-:" + sf[d].sdatam + " was not found in mtz");
      if (!devmcol[d])
 	Bp3Error("Crystal::setwithmtz", "Column for SF- or SI-:" + sf[d].sdevm + " was not found in mtz");
    } 

    if (sf[d].sfmodel.size())
      for (unsigned m          = 0; m < sf[d].sfmodel.size(); m++)
	if (!modelfcol[d][m])
	  Bp3Error("Crystal::setwithmtz", "Column for FC:" + sf[d].sfmodel[m] + " was not found in mtz");

    if (invert)
    {
      string tempda              = sf[d].Getsdatap();
      string tempde              = sf[d].Getsdevp();
      sf[d].Setsdatap(sf[d].Getsdatam());
      sf[d].Setsdevp(sf[d].Getsdevm());
      sf[d].Setsdatam(tempda);
      sf[d].Setsdevm(tempde);
    }
  }  

  // maximum number of reflections in mtz file
  maxref                         = CMtz::MtzNref(MTZIN);
  maxselref                      = maxref;
  resize();

  int ind_xtal, ind_set, ind_col[3];
  if (!CMtz::MtzFindInd(MTZIN,&ind_xtal, &ind_set, ind_col))
    Bp3Error("Crystal::setwithmtz", "HKL could not be found in mtz");

  Setdataset();

  // setup the different cells
  for (unsigned c                = 0; c < cell.size(); c++)
    if (cell[c].a               == 0.0)
    { 	
      unsigned d(dataset[c][0]);
      for (unsigned j            = 0; j < col.size(); j++)
	if (!strcmp(col[j]->label, sf[dataset[c][0]].sdatap.c_str()))
        {
	  xtl                    = CMtz::MtzIxtal(MTZIN,xtalnumber[j]);
	  cell[c].Setcell(xtl->cell);
	}
    }  

  vector<work> highestres(sf.size(), 1000.0);
  vector<work> lowestres(sf.size(), ZERO);
  
  for (unsigned d                = 0; d < sf.size(); d++)
    sf[d].nref                   = 0;

  for (unsigned r                = 0; r < maxref; r++)
  {

    miller[r][0]                 = (int)MTZIN->xtal[ind_xtal]->set[ind_set]->col[ind_col[0]]->ref[r];
    miller[r][1]                 = (int)MTZIN->xtal[ind_xtal]->set[ind_set]->col[ind_col[1]]->ref[r];
    miller[r][2]                 = (int)MTZIN->xtal[ind_xtal]->set[ind_set]->col[ind_col[2]]->ref[r];
    
    // set cell
    for (unsigned c              = 0; c < cell.size(); c++)
      cell[c].Setstolsq(miller,r);

    // set epsilon and centricity
    Setepsilon(r).Setcentric(r);

    if (Sysabs(r))
    {
      if (verbose)
	Bp3Warning("Crystal::setwithmtz","Systematically absent reflection found - it won't be used or phased");
      for (unsigned d            = 0; d < sf.size(); d++)
      {
	sf[d].datap[r]           = NOTUSED;
 	sf[d].devp[r]            = NOTUSED;
	if (sf[d].anomalous)
	{
	  sf[d].datam[r]         = NOTUSED;
	  sf[d].devm[r]          = NOTUSED;
	}
      }
    }
    else
      for (unsigned d            = 0; d < sf.size(); d++)
      {
	if (!CMtz::ccp4_ismnf(MTZIN,col[datapcol[d]]->ref[r]) && (col[devpcol[d]]->ref[r] > ZERO) )
	{
	  // do not store reflections with non-positive sigmas
	  sf[d].datap[r]         = (work)col[datapcol[d]]->ref[r];
	  sf[d].devp[r]          = (work)col[devpcol[d]]->ref[r];
	}
	else
	{
	  sf[d].datap[r]         = NOTUSED;
	  sf[d].devp[r]          = NOTUSED;
	}

	for (unsigned m          = 0; m < sf[d].sfmodel.size(); m++)
	  if (sf[d].sfmodel[m].size())
	    if (!CMtz::ccp4_ismnf(MTZIN,col[modelfcol[d][m]]->ref[r]) &&
		!CMtz::ccp4_ismnf(MTZIN,col[modelpcol[d][m]]->ref[r])   )
	    { 	
	      sf[d].fmodel[r][m] = (work)col[modelfcol[d][m]]->ref[r];
	      sf[d].pmodel[r][m] = (work)col[modelpcol[d][m]]->ref[r];
	    }
	    else
	    {
	      sf[d].fmodel[r][m] = NOTUSED;
	      sf[d].pmodel[r][m] = NOTUSED;
	    }
	
	if (sf[d].Getsfom().size())
	  if (!CMtz::ccp4_ismnf(MTZIN,col[fomcol[d]]->ref[r]))
	    sf[d].fom[r]         = (work)col[fomcol[d]]->ref[r];
	  else
	    sf[d].fom[r]         = NOTUSED;

	if (sf[d].Getsphi().size())
	  if (!CMtz::ccp4_ismnf(MTZIN,col[phicol[d]]->ref[r]))
	    sf[d].phib[r]        = (work)col[phicol[d]]->ref[r];
	  else
	    sf[d].phib[r]        = NOTUSED;

	if (sf[d].Getshla().size())
	  if (!CMtz::ccp4_ismnf(MTZIN,col[hlacol[d]]->ref[r]))
	    sf[d].hla[r]         = (work)col[hlacol[d]]->ref[r];
	  else
	    sf[d].hla[r]         = ZERO;

	if (sf[d].Getshlb().size())
	  if (!CMtz::ccp4_ismnf(MTZIN,col[hlbcol[d]]->ref[r]))
	    sf[d].hlb[r]         = (work)col[hlbcol[d]]->ref[r];
	  else
	    sf[d].hlb[r]         = ZERO;

	if (sf[d].Getshlc().size())
	  if (!CMtz::ccp4_ismnf(MTZIN,col[hlccol[d]]->ref[r]))
	    sf[d].hlc[r]         = (work)col[hlccol[d]]->ref[r];
	  else
	    sf[d].hlc[r]         = ZERO;

	if (sf[d].Getshld().size())
	  if (!CMtz::ccp4_ismnf(MTZIN,col[hldcol[d]]->ref[r]))
	    sf[d].hld[r]         = (work)col[hldcol[d]]->ref[r];
	  else
	    sf[d].hld[r]         = ZERO;

	if (sf[d].anomalous)
	  if (!CMtz::ccp4_ismnf(MTZIN,col[datamcol[d]]->ref[r]) && (col[devmcol[d]]->ref[r] > ZERO) ) 
	  {
	    sf[d].datam[r]       = (work)col[datamcol[d]]->ref[r];
	    sf[d].devm[r]        = (work)col[devmcol[d]]->ref[r];
	  }
	  else if (centric[r] && (sf[d].datap[r] >= NOTUSED) )  
	  {
	    // ***NSP - default in many programs is to have F- == absent for centrics
	    sf[d].datam[r]       = sf[d].datap[r];
	    sf[d].devm[r]        = sf[d].devp[r];
	  }
	  else
	  {
	    sf[d].datam[r]       = NOTUSED;
	    sf[d].devm[r]        = NOTUSED;
	  }
      }
    
    // sigma cut off's and resolution limits

    for (unsigned d              = 0; d < sf.size(); d++)
    {
      if (sf[d].datap[r]         > NOTUSED)
      {
	if (sf[d].datap[r]/sf[d].devp[r] < sf[d].sigmacut)
	  sf[d].datap[r]         = sf[d].devp[r] = NOTUSED;
	if ( (d                  > 0) && sf[0].use(r) )
	  if (fabs(sf[d].datamean(r) - sf[0].datamean(r))/
	      sqrt(sf[d].devmean(r)*sf[d].devmean(r) +
		   sf[0].devmean(r)*sf[0].devmean(r)) < sf[d].isosigmacut)
	  {
	    sf[d].datap[r]       = sf[d].devp[r] = NOTUSED;
	    if (sf[d].anomalous)
	      sf[d].datam[r]     = sf[d].devm[r] = NOTUSED;
	  }
	
      	/*
	if (sf[d].anomalous)
	{
	  if (sf[d].datam[r] > NOTUSED)
	    if (sf[d].datam[r]/sf[d].devm[r] < sf[d].sigmacut)
	      sf[d].datam[r] = sf[d].devm[r] = NOTUSED;
	  if (sf[d].anouse(r) )
	    if (fabs(sf[d].dano(r))/sf[d].sigdano(r) < sf[d].anosigmacut)
	      sf[d].datam[r] = sf[d].devm[r] = sf[d].datap[r] = sf[d].devp[r] = NOTUSED;
	}
	*/
      }
      work resol                 = Getres(d,r);
      // calculate the highest and lowest resolution limits
      if (sf[d].use(r))
      {
 	highestres[d]            = std::min(resol, highestres[d]);
	lowestres[d]             = std::max(resol, lowestres[d]);
      }

      // if the user has specified resolution limits,
      // remove any reflections outside of the boundary 
      if (sf[d].hires           != sf[d].MAXRES)
	if (resol                < sf[d].hires)
	{
	  sf[d].datap[r]         = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]       = NOTUSED;
	}
      if (sf[d].lowres          != sf[d].MINRES)
	if (resol                > sf[d].lowres)
	{
	  sf[d].datap[r]         = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]       = NOTUSED;
	}
      if (sf[d].use(r))
	sf[d].nref++;
    }
    
  }

  if (MTZIN)
    CMtz::MtzFree(MTZIN);  

  for (unsigned d                = 0; d < sf.size(); d++)
  {
    // user hasn't set resolution limits, so set it here.
    if (sf[d].hires             == sf[d].MAXRES)
      sf[d].hires                = highestres[d];
    if (sf[d].lowres            == sf[d].MINRES)
      sf[d].lowres               = lowestres[d];

    if (sf[d].nwave              > 0)
      mad                        = true;
  
    // set nbins, if not inputted by user
    if (!sf[d].nbins)        
      sf[d].Setnbins();
    sf[d].resizebin();
    sf[d].Setshell();
  }

  if (mad)     // set nbins == 5 (if user hasn't set it)
    Setonlyacentrics(true);  // ***NSP

  Setselected(type);

  for (unsigned d                = 0; d < sf.size(); d++)
  {
    sf[d].checkparameters();
    if (clipper)
    {
      //clipper::Message::set_message_level(1);
      double a(cell[sf[d].nxtal].a), b(cell[sf[d].nxtal].b);
      double c(cell[sf[d].nxtal].c), alpha(cell[sf[d].nxtal].alpha);
      double beta(cell[sf[d].nxtal].beta), gamma(cell[sf[d].nxtal].gamma);

      if (spg_trial             == -1)
	spg_trial                = sg.Getnumber();

      clipper::Spacegroup sg_hkl_data(clipper::Spgr_descr(sg.Getnumber()));
      clipper::Spacegroup sg_hkl_data_trial(clipper::Spgr_descr(spg_trial));
      clipper::Cell cell_hkl_data(clipper::Cell_descr(a, b, c, alpha, beta, gamma));
      clipper::Resolution reso_hkl_data(sf[d].hires-0.01);

      sf[d].hkl_list_data_mtz.init(sg_hkl_data, cell_hkl_data, reso_hkl_data, true);
      sf[d].hkl_list_data.init(clipper::Spacegroup(clipper::Spgr_descr(spg_trial)), cell_hkl_data, reso_hkl_data, true);
      //sf[ndata].hkl_list.init(model_base.spacegroup(), model_base.cell(), reso_hkl_data, true);
      sf[d].hkl_list.init(clipper::Spacegroup::p1(), cell_hkl_data, reso_hkl_data, true);

      sf[d].Fna.init(sf[d].hkl_list, sf[d].hkl_list.cell());
      sf[d].Fobs_P1.init(sf[d].hkl_list, sf[d].hkl_list.cell());
      sf[d].Fobs_mtz.init(sf[d].hkl_list_data_mtz, sf[d].hkl_list.cell());
      sf[d].Fobs.init(sf[d].hkl_list_data, sf[d].hkl_list.cell());
      if (phasematch)
	sf[d].Fref.init(sf[d].hkl_list_data, sf[d].hkl_list_data.cell());

      for (unsigned r            = 0; r < maxselref; r++)
	if (selected[r])
	{
	  clipper::HKL_info::HKL_reference_coord rx(sf[d].hkl_list_data_mtz, clipper::HKL(miller[r][0], miller[r][1], miller[r][2]));
	  sf[d].Fobs_mtz[rx.index()].f() = sf[d].datap[r];
	  sf[d].Fobs_mtz[rx.index()].sigf() = sf[d].devp[r];
	}

      model_result.init(clipper::Spacegroup(clipper::Spgr_descr(spg_trial)), cell_hkl_data);
      vector <clipper::Atom> model_base_atom_list(model_base.model().atom_list());
      vector <clipper::Atom> ref_model_atom_list(ref_model.model().atom_list());

      model_base_atom_list.erase(remove_if(model_base_atom_list.begin(), model_base_atom_list.end(), filter_atom), model_base_atom_list.end());
      ref_model_atom_list.erase(remove_if(ref_model_atom_list.begin(), ref_model_atom_list.end(), filter_atom), ref_model_atom_list.end());
      clipper::SFcalc_iso_fft<float>(sf[d].Fna, model_base_atom_list);
	  //for (clipper::HKL_info::HKL_reference_index ih = sf[d].hkl_list.first(); !ih.last(); ih.next() )
	  //	std::cout << std::setw(5) << ih.index() << std::setw(5) << ih.hkl().h() <<   std::setw(5)<< ih.hkl().k() <<  std::setw(5)<< ih.hkl().l() << std::setw(15) << sf[d].Fna[ih].f() << std::setw(15) << sf[d].Fna[ih].phi() << std::endl; 
      if (phasematch)
	clipper::SFcalc_iso_fft<float>(sf[d].Fref, ref_model_atom_list);
    //      sf[d].resizeclipper();
    //
    //  Change dataset to different spacegroup
      for (clipper::HKL_info::HKL_reference_index ih = sf[d].hkl_list.first(); !ih.last(); ih.next() )
	if (!(sf[d].Fobs_mtz[ih.hkl()].missing()))
	  sf[d].Fobs_P1[ih]      = sf[d].Fobs_mtz[ih.hkl()];


      for (clipper::HKL_info::HKL_reference_index ih = sf[d].hkl_list_data.first(); !ih.last(); ih.next() )
	if (!(sf[d].Fobs_P1[ih.hkl()].missing())) {
	  sf[d].Fobs[ih] = sf[d].Fobs_P1[ih.hkl()];
				//std::cout << std::setw(5) << ih.index() << std::setw(5) << ih.hkl().h() <<   std::setw(5)<< ih.hkl().k() <<  std::setw(5)<< ih.hkl().l() << std::setw(15) << sf[d].Fobs[ih].f() << std::setw(15) << sf[d].Fobs[ih].sigf() << std::endl; 
	}	
    }
  }

  if (clipper)
  {
    Setnormalization();
    for (unsigned d              = 0; d < sf.size(); d++)
      sf[d].resizeclipper();
  }

  return *this;
}

Crystal &Crystal::setwithpdb(char *filename)
{
  clipper::MMDBfile mfile;
  clipper::MiniMol model_base_tmp;
  mfile.read_file( filename );
  mfile.import_minimol( model_base_tmp );

  model_base.init(model_base_tmp.spacegroup(), model_base_tmp.cell());
  prepare_pdb_model (model_base, model_base_tmp);

//  for ( int p = 0; p < model_base.size(); p++ )
//   for ( int m = 0; m < model_base[p].size(); m++ )
//    for ( int a = 0; a < model_base[p][m].size(); a++ )
//     std::cout << model_base[p].id()+"\t"+model_base[p][m].id()+"\t"+model_base[p][m].type()+"\t"+model_base[p][m][a].id() << 
//		"\t" << model_base[p][m][a].coord_orth().x() << "\t" << model_base[p][m][a].coord_orth().y() << "\t" << model_base[p][m][a].coord_orth().z() <<
//		"\t" << model_base[p][m][a].occupancy() << "\t" << model_base[p][m][a].u_iso()  << std::endl;

  phasematch = false;
  return *this;
}

Crystal &Crystal::setwithpdb_ref(char *filename)
{
  clipper::MMDBfile mfile;
  mfile.read_file( filename );
  mfile.import_minimol( ref_model );
  phasematch = true;
  return *this;
}

void Crystal::resize()
{
  // allocate memory for vectors holding X-ray data
  
  miller.resize(maxref);

  selected.resize(maxref, true);

  // allocate memory for vectors holding various 
  // information about the particular reflection

  // centricity and phase restriction
  centric.resize(maxref, false);
  centricphase.resize(maxref, ZERO);
  // multiplicity
  epsilon.resize(maxref, 0);
  
  for (unsigned r     = 0; r < maxref; r++)
    miller[r].resize(3);  

  for (unsigned d     = 0; d < sf.size(); d++)
    if (sf[d].anomalous || sf[d].sdatam.size())
      sf[d].resize(maxref, true);
    else
      sf[d].resize(maxref, false);
}

Crystal &Crystal::Setselected(const string type)
{
  // Select the reflections we can refine  
  // and accumulate some stats
  
  unsigned nsel(0), refcount(1), nrej(0);
  if (type                   == "REFINE")
    refcount                  = 2;

  for (unsigned d             = 0; d < sf.size(); d++)
    for (unsigned s           = 0; s < sf[d].nbins; s++)
    {
      sf[d].nshl[s]           = 0;
      sf[d].astolsq[s]        = ZERO;
      sf[d].fovers[s]         = ZERO;
      if (sf[d].anomalous)
      {
	sf[d].anonshl[s]      = 0;
	sf[d].danoovers[s]    = ZERO;
      }
    }    

  for (unsigned r             = 0; r < maxref; r++)
  {
    selected[r]               = true;

    for (unsigned d           = 0; d < sf.size(); d++)
      if (sf[d].datap[r]     >= NOTUSED)
      {
	if (sf[d].datap[r]/sf[d].devp[r] < sf[d].sigmacut)
	  sf[d].datap[r]      = sf[d].devp[r] = NOTUSED;
	if ( (d               > 0) && sf[0].use(r) )
	  if (fabs(sf[d].datamean(r) - sf[0].datamean(r))/
	      sqrt(sf[d].devmean(r)*sf[d].devmean(r) +
		   sf[0].devmean(r)*sf[0].devmean(r)) < sf[d].isosigmacut)
	  {
	    sf[d].datap[r]    = sf[d].devp[r] = NOTUSED;
	    if (sf[d].anomalous)
	      sf[d].datam[r]  = sf[d].devm[r] = NOTUSED;
	    nrej++;
	  }

  	if (sf[d].anomalous)
	{
	  if (sf[d].datam[r] >= NOTUSED)
	    if (sf[d].datam[r]/sf[d].devm[r] < sf[d].sigmacut)
	    {
	      nrej++;
	      sf[d].datam[r]  = sf[d].devm[r] = NOTUSED;
	    }
 	  if (sf[d].anouse(r) )
	  {
	    if (fabs(sf[d].dano(r))/sf[d].sigdano(r) < sf[d].anosigmacut)
	    {
	      nrej++;
	      sf[d].datam[r]  = sf[d].devm[r] = sf[d].datap[r] = sf[d].devp[r] = NOTUSED;
	    }
	    // ***NSP
	    if (fabs(sf[d].dano(r)) > sf[d].anormscut*std::min(sf[d].datam[r],sf[d].datap[r]))
	    {
	      nrej++;
	      sf[d].datap[r]  = sf[d].datamean(r);
	      sf[d].devp[r]   = sf[d].devmean(r);
  	      sf[d].datam[r]  = sf[d].devm[r] = NOTUSED;
	    }	    
	  }
	}
      }

    unsigned counts(0);

    
    for (unsigned d           = 0; d < sf.size(); d++)
      if (sf[d].use(r))
      {
	counts++;
	if (sf[d].anomalous)
	  if (sf[d].anouse(r))
	    counts++;
      }

    if (counts               == 1)
      for (unsigned d         = 0; d < sf.size(); d++)
	if (sf[d].use(r)   && (type == "PHASE") )
	  counts              = (sf[d].sigmah[0] > ZERO) ? 1 : 0;

    if (!centric[r])       
      selected[r]             = ((counts < refcount) || onlycentrics)  ? false : true;
    else
      selected[r]             = ((counts < refcount) || onlyacentrics) ? false : true;
    
    if (selected[r])
    {
      nsel++;
      maxselref               = r;
      // Reset number of shells in each bin
      Setnshl(r);
    }
  }

  printf("Number of reflections selected: %u\n\n",nsel);
  printf("Number of reflections rejected based on sigma and/or rms cutoffs: %u\n\n",nrej);
  
  for (unsigned d           = 0; d < sf.size(); d++)
    for (unsigned s         = 0; s < sf[d].nbins; s++)
    {
      if (sf[d].nshl[s])
      {
	sf[d].astolsq[s]   /= (double)(sf[d].nshl[s]);
	sf[d].fovers[s]    /= (double)(sf[d].nshl[s]);
      }
      else
	sf[d].astolsq[s]    = (ONE/(FOUR*sf[d].Getlores(s)*sf[d].Getlores(s)) +
			       ONE/(FOUR*sf[d].Gethires(s)*sf[d].Gethires(s)))/TWO;
      if (sf[d].anomalous && sf[d].anonshl[s])
	sf[d].danoovers[s] /= (double)(sf[d].anonshl[s]);
    }
  
  return *this;  
}

Crystal &Crystal::Setepsilon(const unsigned r)
{
  // Calculate epsilon factor for a reflection.

  epsilon[r]            = 0;
  
  for (unsigned isym    = 0; isym < sg.NSYM; isym++)
  {
    int hhh             = (sg.symrot[isym][0][0]*miller[r][0] +
                           sg.symrot[isym][0][1]*miller[r][1] +
                           sg.symrot[isym][0][2]*miller[r][2]);
    
    int kkk             = (sg.symrot[isym][1][0]*miller[r][0] +
                           sg.symrot[isym][1][1]*miller[r][1] +
                           sg.symrot[isym][1][2]*miller[r][2]);

    int lll             = (sg.symrot[isym][2][0]*miller[r][0] +
                           sg.symrot[isym][2][1]*miller[r][1] +
                           sg.symrot[isym][2][2]*miller[r][2]);

    if ( (miller[r][0] == hhh) && (miller[r][1] == kkk) && (miller[r][2] == lll) )
      epsilon[r]++;
  }

  if ((sg.lattice == 'A') || (sg.lattice == 'B') ||
      (sg.lattice == 'C') || (sg.lattice == 'I') )
    epsilon[r]         *= 2;
  else if ((sg.lattice == 'R') || (sg.lattice == 'H'))
    epsilon[r]         *= 3;
  else if (sg.lattice  == 'F')
    epsilon[r]         *= 4;

  return *this;
}

Crystal &Crystal::Setcentric(const unsigned r)
{
  // Determines whether a reflection is centric or not.                            
  // Also determines the phase restriction and passes it in phase.                 
  // (ie. phase restriction is phase and phase + pi             

  centricphase[r]       = ZERO;
  centric[r]            = false;
    
  for (unsigned isym    = 0; (isym < sg.NSYMP) && !centric[r]; isym++)
  {
    int hhh             = (sg.symrot[isym][0][0]*miller[r][0] +
                           sg.symrot[isym][0][1]*miller[r][1] +
                           sg.symrot[isym][0][2]*miller[r][2]);
    
    int kkk             = (sg.symrot[isym][1][0]*miller[r][0] +
                           sg.symrot[isym][1][1]*miller[r][1] +
                           sg.symrot[isym][1][2]*miller[r][2]);

    int lll             = (sg.symrot[isym][2][0]*miller[r][0] +
                           sg.symrot[isym][2][1]*miller[r][1] +
                           sg.symrot[isym][2][2]*miller[r][2]);
    
    if ( (miller[r][0] == -hhh) && (miller[r][1] == -kkk) && (miller[r][2] == -lll) )
    {
      centric[r]        = true;
      work rr           = (miller[r][0]*sg.symtran[isym][0] +
			   miller[r][1]*sg.symtran[isym][1] +
			   miller[r][2]*sg.symtran[isym][2]);
      work ss           = (work)( (int)rr );
      centricphase[r]   = (rr >= ss) ?  PI*(rr - ss) : PI*(rr - ss + ONE);  
    }
  }
  return *this;
}

bool Crystal::Sysabs(const unsigned r) const
{
  bool abs(false);
  for (unsigned isym             = 0; (isym < sg.NSYM) && !abs; isym++)
  {
    int hhh                      = (sg.symrot[isym][0][0]*miller[r][0] +
				    sg.symrot[isym][0][1]*miller[r][1] +
				    sg.symrot[isym][0][2]*miller[r][2]);
    
    int kkk                      = (sg.symrot[isym][1][0]*miller[r][0] +
				    sg.symrot[isym][1][1]*miller[r][1] +
				    sg.symrot[isym][1][2]*miller[r][2]);

    int lll                      = (sg.symrot[isym][2][0]*miller[r][0] +
				    sg.symrot[isym][2][1]*miller[r][1] +
				    sg.symrot[isym][2][2]*miller[r][2]);

    if ( (miller[r][0]          == hhh) && (miller[r][1] == kkk) && (miller[r][2] == lll) )
      if (norm(sg.symtran[isym]) > ZERO)
	abs                      = (cos(TWOPI*(miller[r][0]*sg.symtran[isym][0] +
					       miller[r][1]*sg.symtran[isym][1] +
					       miller[r][2]*sg.symtran[isym][2])) < 0.9);
	
  }
  return abs;
}

Crystal &Crystal::Setreso()
{
  vector<work> highestres(sf.size(), 1000.0);
  vector<work> lowestres(sf.size(), ZERO);

  for (unsigned d            = 0; d < sf.size(); d++)
  {
    for (unsigned r          = 0; r < maxref; r++)
    {
      work resol             = Getres(d,r);
      // calculate the highest and lowest resolution limits
      if (sf[d].use(r))
      {
	highestres[d]        = std::min(resol, highestres[d]);
	lowestres[d]         = std::max(resol, lowestres[d]);
      }

      // if the user has specified resolution limits,
      // remove any reflections outside of the boundary 
      if (sf[d].hires       != sf[d].MAXRES)
	if (resol            < sf[d].hires)
	{
	  sf[d].datap[r]     = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]   = NOTUSED;
	}
      if (sf[d].lowres      != sf[d].MINRES)
	if (resol            > sf[d].lowres)
	{
	  sf[d].datap[r]     = NOTUSED;
	  if (sf[d].anomalous)
	    sf[d].datam[r]   = NOTUSED;
	}
    }
    // user hasn't set resolution limits, so set it here.
    if (sf[d].hires         == sf[d].MAXRES)
      sf[d].hires            = highestres[d];
    if (sf[d].lowres        == sf[d].MINRES)
      sf[d].lowres           = lowestres[d];

  }
  return *this;
}

Crystal &Crystal::Setdataset()
{
  // temporarily set the dataset to maximum possible

  dataset.resize(sf.size());
  
  for (unsigned d                     = 0; d < sf.size(); d++)
    dataset[d].resize(sf.size());
  
  unsigned maxxtal(0);
  unsigned maxwave(0);
  
  for (unsigned d                     = 0; d < sf.size(); d++)
  {
    maxxtal                           = std::max(maxxtal, sf[d].nxtal);
    maxwave                           = std::max(maxwave, sf[d].nwave);
    dataset[sf[d].nxtal][sf[d].nwave] = d;
  }
  dataset.resize(maxxtal+1);
  for (unsigned x                     = 0; x < dataset.size(); x++)
    dataset[x].resize(maxwave+1);
  return *this;
}

double Crystal::Calculatebfactor(const unsigned nucleotides, const unsigned residues,
				 const unsigned monomers, const unsigned cbd) 
{

  const double ncarbon((residues*CARBONRES + nucleotides*CARBONNUC)*monomers);
  const double nnitrogen((residues*NITRORES + nucleotides*NITRONUC)*monomers);
  const double noxygen((residues*OXYGENRES + nucleotides*OXYGENNUC)*monomers);
  const double nsulfur(residues*SULFURRES*monomers);
  const double nphosphorus(nucleotides*PHOSNUC*monomers);  
  
  // atomic form-factors: rows: a, b, and f',f",c

  const double cff []         = { 2.310000, 1.020000, 1.588600, 0.865000,
				  20.843899, 10.207500, 0.568700, 51.651199,
				  0.017000, 0.009000, 0.215600};
  
  const double nff []         = {12.212600, 3.132200, 2.012500, 1.166300,
				 0.005700, 9.893300, 28.997499, 0.582600,
				 0.029000, 0.018000, -11.528999};

  const double off []         = { 3.048500, 2.286800, 1.546300, 0.867000,
				  13.277100, 5.701100, 0.323900, 32.908897,
				  0.047000, 0.032000, 0.250800};

  const double pff []         = { 6.434500, 4.179100, 1.780000, 1.490800,
				  1.906700, 27.157000, 0.526000, 68.164497,
				  0.283000, 0.434000, 1.114900};

  const double sff []         = { 6.905300, 5.203400, 1.437900, 1.586300,
				  1.467900, 22.215099, 0.253600, 56.172001,
				  0.319000, 0.557000, 0.866900};


  // calculate absolute scale
  vector<double> sumx(sf.size(), ZERO), sumx2(sf.size(), ZERO);
  vector<double> sumxy(sf.size(), ZERO), sumy(sf.size(), ZERO), sumtot(sf.size(), ZERO);
  
  for (unsigned r              = 0; r < maxref; r++)
    if (selected[r]   && !centric[r])
      for (unsigned d          = 0; d < sf.size(); d++)
      {
	double stolsq(cell[sf[d].nxtal].stolsq[r]);
    
	if (ONE/(sqrt(stolsq)) < 8.0)
	{
	  double scatc(cff[0]*exp(-cff[4]*stolsq) +
		       cff[1]*exp(-cff[5]*stolsq) +
		       cff[2]*exp(-cff[6]*stolsq) +
		       cff[3]*exp(-cff[7]*stolsq) +
		       cff[10]);

	  double scatn(nff[0]*exp(-nff[4]*stolsq) +
		       nff[1]*exp(-nff[5]*stolsq) +
		       nff[2]*exp(-nff[6]*stolsq) +
		       nff[3]*exp(-nff[7]*stolsq) +
		       nff[10]);

	  double scato(off[0]*exp(-off[4]*stolsq) +
		       off[1]*exp(-off[5]*stolsq) +
		       off[2]*exp(-off[6]*stolsq) +
		       off[3]*exp(-off[7]*stolsq) +
		       off[10]);

	  double scats(sff[0]*exp(-sff[4]*stolsq) +
		       sff[1]*exp(-sff[5]*stolsq) +
		       sff[2]*exp(-sff[6]*stolsq) +
		       sff[3]*exp(-sff[7]*stolsq) +
		        sff[10]);

	  double scatp(pff[0]*exp(-pff[4]*stolsq) +
		       pff[1]*exp(-pff[5]*stolsq) +
		       pff[2]*exp(-pff[6]*stolsq) +
		       pff[3]*exp(-pff[7]*stolsq) +
		       pff[10]);

	  double sigmannot((scatc*scatc*ncarbon + scatn*scatn*nnitrogen +
			    scato*scato*noxygen + scats*scats*nsulfur +
			    scatp*scatp*nphosphorus)*sg.NSYMP);

	  unsigned sa(bin(d,r)), sb(sa);
	  double wa(ONE), wb(ZERO);
	  binweights(d, r, sa, sb, wa, wb);
	  double sigman         = wa*sf[d].sigman[sa] + wb*sf[d].sigman[sb];

	  double y(sigman/sigmannot);
	  if (y                 > EPSILON)
	  {
	    y                   = -log(y);
	    sumx[d]            += TWO*stolsq;
	    sumx2[d]           += FOUR*stolsq*stolsq;
	    sumxy[d]           += TWO*y*stolsq;
	    sumy[d]            += y;
	    sumtot[d]          += ONE;
	  }
	}
      }

  vector<double> kscale(sf.size(), ZERO), wilsonb(sf.size(), ZERO);
  
  for (unsigned d             = 0; d < sf.size(); d++)
  {
    printf("Dataset %s:\n",sf[d].name.c_str());

    double denom              = sumtot[d]*sumx2[d] - sumx[d]*sumx[d];
    kscale[d]                 = exp((sumy[d]*sumx2[d] - sumx[d]*sumxy[d])/denom);
    wilsonb[d]                = ((sumtot[d]*sumxy[d] - sumx[d]*sumy[d])/denom);

    if (verbose              || rscale)
    {
      if (sf[d].type         == "INTENSITY")
	printf("Absolute scale factor from Wilson plot (Intensities) %f\n", kscale[d]);
      else
      {
	kscale[d]             = sqrt(kscale[d]);
	printf("Absolute scale factor from Wilson plot (Amplitudes) %f\n", kscale[d]);
      }
    }
    wilsonb[d]                = std::min(100.0,std::max(wilsonb[d],5.0));

    printf("B-factor from Wilson plot is %7.3f\n\n", wilsonb[d]);
  }
  
  // likelihood bfactor 
  vector<double> bfactor(sf.size(), ZERO), sumb(sf.size(), ZERO);

  for (unsigned r             = 0; r < maxref; r++)
    for (unsigned d           = 0; d < sf.size(); d++)
    {
      double stolsq(cell[sf[d].nxtal].stolsq[r]);
    
      if (selected[r]   && !centric[r] && (ONE/(sqrt(stolsq)) < 8.0) && sf[d].use(r))
      {
	double scatc(cff[0]*exp(-cff[4]*stolsq) +
		     cff[1]*exp(-cff[5]*stolsq) +
		     cff[2]*exp(-cff[6]*stolsq) +
		     cff[3]*exp(-cff[7]*stolsq) +
		     cff[10]);

	double scatn(nff[0]*exp(-nff[4]*stolsq) +
		     nff[1]*exp(-nff[5]*stolsq) +
		     nff[2]*exp(-nff[6]*stolsq) +
		     nff[3]*exp(-nff[7]*stolsq) +
		     nff[10]);

	double scato(off[0]*exp(-off[4]*stolsq) +
		     off[1]*exp(-off[5]*stolsq) +
		     off[2]*exp(-off[6]*stolsq) +
		     off[3]*exp(-off[7]*stolsq) +
		     off[10]);

	double scats(sff[0]*exp(-sff[4]*stolsq) +
		     sff[1]*exp(-sff[5]*stolsq) +
		     sff[2]*exp(-sff[6]*stolsq) +
		     sff[3]*exp(-sff[7]*stolsq) +
		     sff[10]);

	double scatp(pff[0]*exp(-pff[4]*stolsq) +
		     pff[1]*exp(-pff[5]*stolsq) +
		     pff[2]*exp(-pff[6]*stolsq) +
		     pff[3]*exp(-pff[7]*stolsq) +
		     pff[10]);

	double sigmannot((scatc*scatc*ncarbon + scatn*scatn*nnitrogen +
			  scato*scato*noxygen + scats*scats*nsulfur +
			  scatp*scatp*nphosphorus)*sg.NSYMP);
	
	double y(sf[d].datamean(r)*kscale[d]/(epsilon[r]*sigmannot));

	if (sf[d].type       == "AMPLITUDE")
	  y                  *= kscale[d]*sf[d].datamean(r);

	if ( (y               > DSMALL) && (y < ONE) )
	{
	  bfactor[d]         += -log(y)/(TWO*stolsq);
	  sumb[d]            += ONE;
	}
      }
    }

  for (unsigned d             = 0; d < sf.size(); d++)
  {
    bfactor[d]               /= sumb[d];
    bfactor[d]                = std::min(200.0,std::max(5.0,bfactor[d]));
    printf("Dataset %s:\n",sf[d].name.c_str());
    printf("An estimated likelihood-based overall B factor is %7.3f\n\n",bfactor[d]);
  }  

  double cbfactor(bfactor[cbd]);
       
  if (wilsonb[cbd]           != 5.0)
    if ( (wilsonb[cbd]       < bfactor[cbd]) || (cbfactor == 5.0) )
      cbfactor                = wilsonb[cbd];

  if (rscale)
  {
    for (unsigned r           = 0; r < maxref; r++)
      for (unsigned d         = 0; d < sf.size(); d++)
	if (sf[d].use(r))
	  if (sf[d].anomalous)
	    if (sf[d].anouse(r))
	    {
	      sf[d].datap[r] *= kscale[d];
	      sf[d].devp[r]  *= kscale[d];
	      sf[d].datam[r] *= kscale[d];
	      sf[d].devm[r]  *= kscale[d];
	    }
	    else if (sf[d].usep(r))
	    {
	      sf[d].datap[r] *= kscale[d];
	      sf[d].devp[r]  *= kscale[d];
	    }
	    else
	    {
	      sf[d].datam[r] *= kscale[d];
	      sf[d].devm[r]  *= kscale[d];
	    }
	  else
	  {
	    sf[d].datap[r]   *= kscale[d];
	    sf[d].devp[r]    *= kscale[d];
	  }
  }
  
  return cbfactor;
}

void Crystal::Calculatematthews(const unsigned nucleotides, const unsigned residues, double weight, unsigned *nmon, double *solv) const
{
  vector<double> mat, solvent;
  vector<unsigned> monomer;

  const double carbon((residues*CARBONRES + nucleotides*CARBONNUC)*12.0);
  const double nitro((residues*NITRORES + nucleotides*NITRONUC)*14.0);
  const double oxygen((residues*OXYGENRES + nucleotides*OXYGENNUC)*16.0);
  const double sulfur(residues*SULFURRES*32.0);
  const double phos(nucleotides*PHOSNUC*30.0);  
  const double hydro(residues*HYDRORES + nucleotides*HYDRONUC);  

  if (weight           == ZERO)
    weight              = carbon + nitro + oxygen + sulfur + phos + hydro;
  
  double yint(ZERO), xcint(ZERO), wtint(ZERO), aint(ZERO), sint(ZERO);

  for (unsigned i       = 1; i < 250; i++)
  {
    double matthews     = cell[0].volume/(weight*(double)(sg.NSYM*i));
    double solv         = ONE- 1.23/matthews;
    if (solv            > 0.0) 
    {
      mat.push_back(matthews);
      solvent.push_back(solv);
      monomer.push_back(i);
    }
  }
  
  if (nucleotides      == 0)
  {
    // protein only numbers

    const double res[]  = {1.199,  1.501,  1.650,  1.801,  1.901,  2.001,
			   2.201,  2.401,  2.601,  2.801,  3.101,  3.501,
			   5.001,  5.001,  5.001};
    
    const double yo[ ]  = {0.085,  0.312,  0.400,  0.503,  0.597,  0.729,
			   1.052,  1.781,  2.852,  3.386,  3.841,  4.281,
			   4.592,  1.503,  0.257};
  
    const double xc[ ]  = {2.052,  2.102,  2.122,  2.132,  2.140,  2.155,
			   2.171,  2.182,  2.191,  2.192,  2.205,  2.211,
			   2.210,  2.256,  2.324};

    const double wt[ ]  = {0.213,  0.214,  0.220,  0.225,  0.226,  0.231,
			   0.236,  0.241,  0.242,  0.244,  0.244,  0.244,
			   0.245,  0.446,  0.327};
  
    const double a[ ]   = {28.38,  102.7,  187.5,  339.3,  434.1,  540.5,
			   686.2,  767.9,  835.9,  856.9,  854.0,  849.6,
			   846.7,  136.6,  47.10};

    const double s[ ]   = {0.953,  0.807,  0.775,  0.702,  0.648,  0.640,
			   0.635,  0.589,  0.584,  0.542,  0.500,  0.485,
			   0.480,  1.180,  0.466};

     work hires(1000.0);
     for (unsigned d    = 0; d < sf.size(); d++)
       hires            = std::min(hires,sf[d].hires);
    
     for (unsigned i    = 0; i < 13; i++)
       if (hires        < res[i])
	 if ( (i       == 0) || (i == 12) )
	 {
	   yint         =  yo[i];
	   xcint        =  xc[i];
	   wtint        =  wt[i];
	   aint         =  a[i];
	   sint         =  s[i];
	 }
	 else if (i     < 12)
	 {
	   double c     = res[i] - res[i-1];
	   double wa    = (res[i] - hires)/c;
	   double wb    = (hires - res[i-1])/c;
	   yint         =  wa*yo[i-1] + wb*yo[i];
	   xcint        =  wa*xc[i-1] + wb*xc[i];
	   wtint        =  wa*wt[i-1] + wb*wt[i];
	   aint         =  wa*a[i-1] + wb*a[i];
	   sint         =  wa*s[i-1] + wb*s[i];
	 }
     
  }
  else if (residues    == 0)
  {
    // dna only numbers
    yint                = 1.503;
    xcint               = 2.256;
    wtint               = 0.446;
    aint                = 136.6;
    sint                = 1.180;
  }
  else
  {
    // dna/protein complex numbers
    yint                = 0.257;
    xcint               = 2.324;
    wtint               = 0.327;
    aint                = 47.10;
    sint                = 0.466;
  }

    double maxprob(ZERO);
    double totprob(ZERO);
    vector<double> vmprob(mat.size(), ZERO);
    unsigned maxindex(0);
    for (unsigned i     = 0; i < mat.size(); i++)
    {
      double z          = (mat[i] - xcint)/wtint;
      vmprob[i]         = yint + aint*(exp(-exp(-z)-z*sint+1));
      totprob          += vmprob[i];
      if (vmprob[i]     > maxprob)
      {
	maxprob         = vmprob[i];
	maxindex        = i;
      }
    }
    
    printf(" $TABLE: Matthews coefficient, solvent content, and probability vs number of monomers\n");
    printf("$GRAPHS:");
    printf(" Matthews coefficient vs number of monomers:N:1,2:\n");
    printf("       : Solvent content vs number of monomers:N:1,3:\n");
    printf("       : Kantardjieff and Rupp probability vs number of monomers:N:1,4:\n");
    printf("$$\n Monomers  Matthews-Coeff  Solvent-content    Probability  $$\n$$\n");

    for (unsigned i     = 0; i < mat.size(); i++)
      printf("   %2u        %6.3f           %7.3f          %7.3f\n",monomer[i], mat[i],
	     (ONE-1.23/mat[i]), vmprob[i]/totprob);

    printf("$$\n\n");


    if (*nmon          == 0)
    {
      *nmon             = monomer[maxindex];
      printf("Setting the number of monomers to %2u\n\n", monomer[maxindex]);
      if (*solv       == ZERO)
      {
	*solv          = (ONE-1.23/mat[maxindex]);
	printf("Setting the solvent content to %5.3f\n\n",(ONE-1.23/mat[maxindex]));
      }
      else if ( fabs(*solv - (ONE-1.23/mat[maxindex])) > 0.25*(*solv) )
	Bp3Error("Crystal::Calculatematthews","Input solvent content is 25 % greater than the one calculated here!",!warn);
    }
    else
    {
      unsigned input(mat.size());
      for (unsigned i = 0; i < mat.size(); i++)
	if (*nmon    == monomer[i])
	  input       = i;

      if (input      == mat.size())
	Bp3Error("Crystal::Calculatematthews","Highly improbable number of monomers inputted.",!warn);
      
      if (input      == maxindex)
	printf("Inputted number of monomers correspond to maximal Kantardjieff/Rupp probability\n\n");
      
      if (vmprob[input]/totprob < 0.005)
	Bp3Warning("Crystal::Calculatematthews","Improbable (but still possible) number of monomers inputted.");

      if (*solv       == ZERO)
      {
	*solv          = (ONE-1.23/mat[input]);
	printf("Setting the solvent content to %5.3f.\n\n",(ONE-1.23/mat[input]));
      }
      else if ( fabs(*solv - (ONE-1.23/mat[input])) > 0.25*(*solv) ) 
	Bp3Error("Crystal::Calculatematthews","Input solvent content is 25 % greater than the one calculated here!",!warn);
    }
}

Crystal &Crystal::Setnshl(const unsigned r)
{
  for (unsigned d             = 0; d < sf.size(); d++)
    if (sf[d].use(r))
    {
      unsigned s(bin(d,r));
      sf[d].nshl[s]++;
      sf[d].astolsq[s]       += cell[sf[d].nxtal].stolsq[r];
      sf[d].fovers[s]      += sf[d].foversigma(r);
      if (sf[d].anomalous)
	if (sf[d].anouse(r))
	{
	  sf[d].anonshl[s]++;
	  sf[d].danoovers[s] += fabs(sf[d].dano(r))/sf[d].sigdano(r);
	}
      
    }
  return *this;
}

void Crystal::binweights(const unsigned d, const unsigned r, unsigned &sa,
			 unsigned &sb, double &wa, double &wb) const
{
  double stolsq(cell[sf[d].nxtal].stolsq[r]);
  sa            = between((unsigned)0, (unsigned) (work)((stolsq-sf[d].lowsq)*sf[d].vshell),
			    sf[d].nbins-1);

  if (sa       != 0)
  {
    if (stolsq <= sf[d].astolsq[sa])
      sb        = (sa > 0)             ? sa-1 : sa+1;
    else
      sb        = (sa < sf[d].nbins-1) ? sa+1 : sa-1;
  
    double c    = sf[d].astolsq[sb] - sf[d].astolsq[sa];

    wa          = (sf[d].astolsq[sb] - stolsq)/c;
    wb          = (stolsq - sf[d].astolsq[sa])/c;
  }
}

Crystal &Crystal::Setnormalization()
{
  // Calculate normalization factors.

  vector<vector<double> > sumweight(sf.size()), sumweightref(sf.size());
  vector<vector<double> > sumweightano(sf.size());

  heavyref                         = false;
  if (sf[0].sigmah.size())
    heavyref                       = (sf[0].sigmah[0] > ZERO);
  for (unsigned d                  = 0; d < sf.size()  ; d++)
  {
    if (sf[d].sigmah.size())
      if ( (sf[d].sigmah[0]       <= 0.0) && (nd != -1) )
	nd                         = d;

    sumweight[d].resize(sf[d].nbins, ZERO);
    sumweightref[d].resize(sf[d].nbins, ZERO);
    sumweightano[d].resize(sf[d].nbins, ZERO);
    for (unsigned s                = 0; s < sf[d].nbins; s++)
    {
      sf[d].sigman[s]              = ZERO;
      sf[d].sigmap[s]              = ZERO;
      sf[d].sigmanref[s]           = ZERO;
      sf[d].sigmadano[s]           = ZERO;
    }
  }

  for (unsigned r                  = 0; r < maxselref; r++)
    if (selected[r])
    {
      double weight                = (centric[r]) ? ONE : TWO;
      double eps                   = (double)epsilon[r];

      for (unsigned d              = 0; d < sf.size(); d++)
      {
	unsigned s                 = bin(d,r);

	if (sf[d].use(r) && (sf[d].datamean(r) > ZERO))
	{
	  double data(sf[d].datamean(r));
	  if (sf[d].type          == "AMPLITUDE")
	    data                  *= sf[d].datamean(r);
	    
	  if (clipper)
	  {
            clipper::HKL hkl_val(miller[r][0], miller[r][1], miller[r][2]);
	    sf[d].sigman[s]       += weight*(double)(sf[d].Fobs_mtz[hkl_val].f()*sf[d].Fobs_mtz[hkl_val].f())/eps;
	    sf[d].sigmap[s]       += weight*(double)(sf[d].Fna[hkl_val].f()*sf[d].Fna[hkl_val].f())/eps;
		//std::cout << std::setw(5) << hkl_val.h() <<  std::setw(5) << hkl_val.k() << std::setw(5) << hkl_val.l() << std::setw(10) << sf[d].Fobs_mtz[hkl_val].f() << std::setw(10) << sf[d].Fna[hkl_val].f() <<  std::setw(10) << sf[d].datamean(r) << std::endl;
	    
	  }
	  else
	  {
	    sf[d].sigman[s]       += weight*(double)(data)/eps;
	    if (sf[d].sfmodel.size()) // ***NSP should allow multiple model sigmap
	      if (sf[d].fmodel[r][0] > NOTUSED)
		sf[d].sigmap[s]   += weight*(double)(sf[d].fmodel[r][0]*sf[d].fmodel[r][0])/eps;
	  }
	    
	  sumweight[d][s]         += weight;

	  // ***NSP
	  if (sf[d].anomalous)
	    if (sf[d].anouse(r)   && !centric[r])
	    {
	      sf[d].sigmadano[s]  += weight*(double)(TWO*sf[d].dano(r)*sf[d].dano(r))/eps;
	      sumweightano[d][s]  += weight;
	    }
	}
	if (sf[0].use(r))
	{
	  double data(sf[0].datamean(r));
	  if (sf[d].type          == "AMPLITUDE")
	    data                  *= sf[0].datamean(r);
	  sf[d].sigmanref[s]      += weight*(double)(data)/eps;
	  sumweightref[d][s]      += weight;
	}
      }
    }

  for (unsigned d                  = 0; d < sf.size(); d++)
    for (unsigned s                = 0; s < sf[d].nbins; s++)
    {    
      if (sumweight[d][s]          > ZERO)
      {
	sf[d].sigman[s]           /= sumweight[d][s];
	sf[d].sigmap[s]           /= sumweight[d][s];
      }
    
      if (sumweightref[d][s]       > ZERO)   
	sf[d].sigmanref[s]        /= sumweightref[d][s];
      
      if (sf[d].anomalous)
	if (sumweightano[d][s]     > ZERO)
	  sf[d].sigmadano[s]      /= sumweightano[d][s];
      
    }

  return *this;
}

vector<float> Crystal::Getcell(const unsigned d) const
{
  // return reference cell dimensions
  vector<float> returncell(6);

  returncell[0] = (float) cell[sf[d].nxtal].a;
  returncell[1] = (float) cell[sf[d].nxtal].b;
  returncell[2] = (float) cell[sf[d].nxtal].c;
  returncell[3] = (float) cell[sf[d].nxtal].alpha;
  returncell[4] = (float) cell[sf[d].nxtal].beta;
  returncell[5] = (float) cell[sf[d].nxtal].gamma;
   
  return returncell;
}

Crystal &Crystal::SetLuzzati()
{
  // Get initial guess for Luzzati parameters... if required

  for (unsigned d = 0; d < sf.size(); d++)
    if (sf[d].sfmodel.size())
      Guessdmod();  

  Guessdluz().Guesssdluz();
  
  Guessadluz();

  // .Guesseluz() .Guesssdluzref()
  return *this;
}

Crystal &Crystal::Guessdluz()
{
  // Get initial guess for dluz, if user hasn't supplied values for it.
  // The guess is scaled sqrt of the correlation coefficient between
  // E_n^2 and E_0^2  

  sf[0].noniso                      = false;
  
  for (unsigned d                   = 0; d < sf.size(); d++)
  {
    if (!sf[d].dluz.size())  // set dluz, if user hasn't set it
    {
      sf[d].dluz.resize(sf[d].nbins, ONE);

      if (d                        != 0)
      {
	vector<double> sumfsquared(sf[d].nbins, ZERO);
	vector<double> sumfforth(sf[d].nbins, ZERO);
	vector<double> sumfmixed(sf[d].nbins, ZERO);
	vector<double> sumfsquared0(sf[d].nbins, ZERO);
	vector<double> sumfforth0(sf[d].nbins, ZERO);
	vector<unsigned> NinShl(sf[d].nbins, 0);

	for (unsigned r             = 0; r < maxselref; r++)
	  if (selected[r])
	  {
	    double eps              = (double)(epsilon[r]);
	    if (sf[0].use(r) && sf[d].use(r))
	    {
	      unsigned s            = bin(d,r);
	      NinShl[s]++;
	      if (sf[d].type       == "AMPLITUDE")
	      {
		double eobsref      = sf[0].datamean(r)/sqrt(eps*sf[d].sigmanref[s]);
		double eobsd        = sf[d].datamean(r)/sqrt(eps*sf[d].sigman[s]);
		sumfsquared0[s]    += eobsref*eobsref;
		sumfforth0[s]      += eobsref*eobsref*eobsref*eobsref;
		sumfsquared[s]     += eobsd*eobsd;
		sumfforth[s]       += eobsd*eobsd*eobsd*eobsd;
		sumfmixed[s]       += eobsref*eobsd*eobsref*eobsd;
	      }
	      else if (sf[d].type  == "INTENSITY")
	      {
		double eobsref      = fabs(sf[0].datamean(r))/(eps*sf[d].sigmanref[s]);
		double eobsd        = fabs(sf[d].datamean(r))/(eps*sf[d].sigman[s]);
		sumfsquared0[s]    += eobsref;
		sumfforth0[s]      += eobsref*eobsref;
		sumfsquared[s]     += eobsd;
		sumfforth[s]       += eobsd*eobsd;
		sumfmixed[s]       += eobsref*eobsd;		
	      }
	    }
	  }

	for (unsigned s             = 0; s < sf[d].nbins ; s++)
	{
	  double luzz               = ((double)(NinShl[s])*sumfmixed[s] -
				       sumfsquared[s]*sumfsquared0[s]);
	  double denom              = (((double)(NinShl[s])*sumfforth[s] -
					sumfsquared[s]*sumfsquared[s])*
				       ((double)(NinShl[s])*sumfforth0[s] -
					sumfsquared0[s]*sumfsquared0[s]));
	  luzz                      = (denom > ZERO) ? sqrt(luzz/sqrt(denom)) : 0.01;
	  sf[d].dluz[s]             = ((sf[d].sigmanref[s] > ZERO) ? 
				       luzz*sqrt(sf[d].sigman[s]/
						 sf[d].sigmanref[s]): ONE);
	}

	printf(" $TABLE: Dataset %s Isomorphic correlation vs. Resolution\n",
               sf[d].name.c_str());
        printf("$GRAPHS: Isomorphic correlation vs. Res:N:4,7:\n");
        printf("$GRAPHS: Isomorphic correlation vs. STOL2:N:5,7:\n");
        printf("$$\n Bin   LoRes  HiRes    Res    STOL2   Refls   Isomorphism $$\n$$\n");
	double luzzabovethres(ZERO), totalluzz(ZERO);
	unsigned totalnref(0);
	for (unsigned s             = 0; s < sf[d].nbins; s++)
	  if (sf[d].sigman[s]       > ZERO)
	  {
            double luzz             = sf[d].dluz[s]*sqrt(sf[d].sigmanref[s]/
							 sf[d].sigman[s]);
            printf(" %2d  %6.2f  %6.2f  %6.2f  %7.5f %5i  %10.5f\n",s+1,sf[d].Getlores(s),
                   sf[d].Gethires(s),sf[d].Getavres(s),sf[d].astolsq[s],NinShl[s],luzz);
	    if (luzz                > 0.999)
	      luzzabovethres       += ONE;
	    totalluzz              += luzz*NinShl[s];
	    totalnref              += NinShl[s];
          }
	printf("$$\nTOTAL\n");
	printf("     %6.2f  %6.2f                  %5i  %10.5f",sf[d].Getlores(0),
	       sf[d].Gethires(sf[d].nbins-1), totalnref, totalluzz/totalnref);

        printf("\n\n");

	
	 sf[d].noniso                = ((totalluzz/totalnref < 0.99) || (double(luzzabovethres/sf[d].nbins) < 0.5));

	if (sf[d].noniso)
	  printf("Non-isomorphism errors will be considered for Dataset %s\n\n",sf[d].name.c_str());
	else
	  printf("Non-isomorphism errors will not be considered for Dataset %s\n\n",sf[d].name.c_str());
      }
    }
  }
  return *this;
}

Crystal &Crystal::Guessdmod()
{
  // Get initial guess for dluz, if user hasn't supplied values for it.
  // The guess is scaled sqrt of the correlation coefficient between
  // E_C^2 and E_0^2  

  for (unsigned d                   = 0; d < sf.size(); d++)
  {
    if (!sf[d].dmod.size())  // set dluz, if user hasn't set it
    {
      vector<vector<double> > sumfmixed(sf[d].nbins);
      vector<vector<double> > sumfcsquared(sf[d].nbins), sumfcforth(sf[d].nbins);
      sf[d].dmod.resize(sf[d].nbins);

      for (unsigned s               = 0; s < sf[d].nbins; s++)
      {
	sf[d].dmod[s].resize(sf[d].sfmodel.size(), ONE);      
	sumfmixed[s].resize(sf[d].sfmodel.size(), ZERO);
	sumfcsquared[s].resize(sf[d].sfmodel.size(), ZERO);
	sumfcforth[s].resize(sf[d].sfmodel.size(), ZERO);
      }
      
      vector<unsigned> NinShl(sf[d].nbins, 0);

      vector<double> sumfosquared(sf[d].nbins, ZERO);
      vector<double> sumfoforth(sf[d].nbins, ZERO);

      for (unsigned r               = 0; r < maxselref; r++)
	if (selected[r])
	{
	  double eps                = (double)(epsilon[r]);
	  if (sf[d].use(r))
	  {
 	    unsigned s              = bin(d,r);
	    NinShl[s]++;
	    double eobsref;
 	    if (sf[d].type         != "AMPLITUDE")
	      Bp3Error("Crystal::Guessdmod","Amplitudes are required for Guessdmod");
	    
	    eobsref                 = sf[0].datamean(r)/sqrt(eps*sf[d].sigmanref[s]);
	    sumfosquared[s]        += eobsref*eobsref;
	    sumfoforth[s]          += eobsref*eobsref*eobsref*eobsref;

	      for (unsigned m       = 0; m < sf[d].sfmodel.size(); m++)
	      {
		double ecalcm       = sf[d].fmodel[r][m]/sqrt(eps*sf[d].sigmap[s]);
		sumfcsquared[s][m] += ecalcm*ecalcm;
		sumfcforth[s][m]   += ecalcm*ecalcm*ecalcm*ecalcm;
		sumfmixed[s][m]    += eobsref*ecalcm*eobsref*ecalcm;
	      }
	    }
	  }

	for (unsigned s             = 0; s < sf[d].nbins ; s++)
	  for (unsigned m           = 0; m < sf[d].sfmodel.size(); m++)
	  {
	    double luzz             = std::max((double)(NinShl[s])*sumfmixed[s][m] -
 				               sumfcsquared[s][m]*sumfosquared[s], ZERO);
	    double denom            = (((double)(NinShl[s])*sumfoforth[s] -
					sumfosquared[s]*sumfosquared[s])*
				       ((double)(NinShl[s])*sumfcforth[s][m] -
					sumfcsquared[s][m]*sumfcsquared[s][m]));
 	    luzz                    = (denom > ZERO) ? sqrt(luzz/sqrt(denom)) : 0.01;
	    sf[d].dmod[s][m]        = ((sf[d].sigmap[s] > ZERO) ? 
				       luzz*sqrt(sf[d].sigman[s]/
						 sf[d].sigmap[s]): ONE);
	  }

	for (unsigned m             = 0; m < sf[d].sfmodel.size(); m++)
	{
	  printf(" $TABLE: Dataset %s Obs/Model correlation vs. Resolution\n",
		 sf[d].name.c_str());
	  printf("$GRAPHS: Obs/Model correlation vs. Res:N:4,7:\n");
	  printf("$GRAPHS: Obs/Model correlation vs. STOL2:N:5,7:\n");
	  printf("$$\n Bin   LoRes  HiRes    Res    STOL2   Refls   Isomorphism $$\n$$\n");
	  double luzzabovethres(ZERO), totalluzz(ZERO);
	  unsigned totalnref(0);
	   for (unsigned s          = 0; s < sf[d].nbins; s++)
	     if (sf[d].sigmap[s]    > ZERO)
	     {
	       double luzz          = std::max(0.01, sf[d].dmod[s][m]*
			                       sqrt(sf[d].sigmap[s]/sf[d].sigman[s]));
	       printf(" %2d  %6.2f  %6.2f  %6.2f  %7.5f %5i  %10.5f\n",s+1,sf[d].Getlores(s),
		      sf[d].Gethires(s),sf[d].Getavres(s),sf[d].astolsq[s],NinShl[s],luzz);
	       totalluzz           += luzz*NinShl[s];
	       totalnref           += NinShl[s];
	     }
	   printf("$$\nTOTAL\n");
	   printf("     %6.2f  %6.2f                  %5i  %10.5f",sf[d].Getlores(0),
		  sf[d].Gethires(sf[d].nbins-1), totalnref, totalluzz/totalnref);
	   printf("\n\n");
	   printf(" Overall Correlation is %.3f\n\n",totalluzz/totalnref);
	}
    }
  }
  return *this;
}

Crystal &Crystal::Guesseluz()
{
  // initial guess for eluz - if not supplied by the user.

  for (unsigned d = 0; d < sf.size(); d++)
    if (!sf[d].eluz.size())
      sf[d].eluz.resize(sf[d].nbins, ONE);

  return *this;
}

Crystal &Crystal::Guesssdluz()
{
  for (unsigned d = 0; d < sf.size(); d++)
    if (!sf[d].sdluz.size() && sf[d].anomalous)
      sf[d].sdluz.resize(sf[d].nbins, sf[d].defsdluz);  
  
  return *this;
}

/*
Crystal &Crystal::Guesssdluzref()
{
  for (unsigned d = 0; d < sf.size(); d++)
    if (!sf[d].sdluzref.size())
      sf[d].sdluzref.resize(sf[d].nbins, sf[d].defsdluz);
  return *this;
}
*/

Crystal &Crystal::Guessadluz()
{
  // initial guess for adluz - if not supplied by the user.

  for (unsigned d = 0; d < sf.size(); d++)
    if (!sf[d].adluz.size() && sf[d].anomalous) 
      sf[d].adluz.resize(sf[d].nbins, 0.85);

  return *this;
}

void Crystal::printscale() const
{
  printf("Scale Parameters\n");

  for (unsigned d     = 0; d  < sf.size() ; d++)
  {
    printf("Dataset: %s\n", sf[d].Getname().c_str());
    printf("K = %10.5f   B = %7.3f\n\n",sf[d].kscale,sf[d].biso);
  }
  
  printf("\n");
  fflush(stdout);
}

void Crystal::printLuzzati(const string targ, const int cyc) const
{
  if (targ           == "MSRS" || targ == "PSAD" || targ == "PSD2" )
  {
    printf("$$\n Bin   HiRes  LoRes   STOL2");
    for (unsigned p   = 0; p < pluz.size(); p++)
      printf("        PLUZ%u",p+1);        
    printf(" $$\n$$\n");
    unsigned d(0);
    if (targ         == "MSRS")
      d++;
    for (unsigned s   = 0; s < sf[d].nbins; s++)
    {
      printf(" %2d  %6.2f  %6.2f  %7.5f",s+1,sf[d].Getlores(s),sf[d].Gethires(s),sf[d].astolsq[s]);
      for (unsigned p = 0; p < pluz.size(); p++)
        printf("   %10.5f",pluz[p][s]);
      printf("\n");
    }
    printf("$$\n\n");
    fflush(stdout);
  }
  else
    for (unsigned d   = 0; d  < sf.size() ; d++)
      if ( ((d       == 0) && (sf[0].anomalous)) || (d != 0) || (sf[d].sfmodel.size()) ) 
	sf[d].printLuzzati(targ, cyc,verbose); 
}

void Crystal::print() const
{
  printf("Crystal information\n");
  printf("Maximum number of reflections: %u\n",maxref); 
  printf("Number of datasets: %u\n\n",(unsigned)(sf.size()));

  sg.print(verbose);

  for (unsigned d   = 0; d < sf.size(); d++)
  {
    printf("Dataset %2u\n",d+1);
    if (sf[d].anomalous)
      printf("Anomalous data present\n");    
    fflush(stdout);
    cell[sf[d].nxtal].print(verbose);
    sf[d].printsignal();
    sf[d].printnorm();
  }
  // printLuzzati();
}
