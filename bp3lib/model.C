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
#include <stdlib.h>
#include "model.h"

// Scatter methods

Scatter::Scatter(const string namein, const unsigned xtalin, const unsigned nwavein)
{
  xtal           = xtalin;
  nwave          = nwavein;
  // resize vectors
  resize();

  Setname(namein);

  // if we have a name, set the normal scattering factors

  if (name         != "")
    Setabc();
  
  // if there is just one wavelength, the default is CuKa 
  // radiation (which is set in Setabc) and not to refine 
  // fp or fpp

  if (nwave        == 1)
  {
    refinefp[0]     = false;
    refinefpp[0]    = false;
  }
  else
  {
    for (unsigned w = 0; w < nwave; w++)
    {
      refinefp[w]   = true;
      refinefpp[w]  = true;
    }
  }
}

Scatter &Scatter::resize()
{
  fp.resize(nwave);
  fpp.resize(nwave);
  refinefp.resize(nwave);
  refinefpp.resize(nwave);

  return *this;
}

Scatter &Scatter::Setfp(const unsigned w, const double value)
{
  if (w   < nwave)
    fp[w] = value;
  else
    Bp3Error("Scatter::Setfp", "Size of input int too large");

  return *this;
}

Scatter &Scatter::Setfp(const vector<double> fpin)
{
  if (nwave == fpin.size())
    fp       = fpin;
  else
    Bp3Error("Scatter::Setfp", "wrong size of input vector");
  
  return *this;
}

Scatter &Scatter::Setfpp(const unsigned w, const double value)
{
  if (w    < nwave)
    fpp[w] = value;
  else
    Bp3Error("Scatter::Setfpp", "Size of input int too large");

  return *this;
}

Scatter &Scatter::Setfpp(const vector<double> fppin)
{
  if (nwave == fppin.size())
    fpp      = fppin;
  else
    Bp3Error("Scatter::Setfpp", "wrong size of input vector");
  
  return *this;
}

Scatter &Scatter::Setrefinefp(const unsigned w, const bool value)
{
  if (w         < nwave)
    refinefp[w] = value;
  else
    Bp3Error("Scatter::Setrefinefp", "Size of input int too large");

  return *this;
}

Scatter &Scatter::Setrefinefpp(const unsigned w, const bool value)
{
  if (w          < nwave)
    refinefpp[w] = value;
  else
    Bp3Error("Scatter::Setrefinefpp", "Size of input int too large");

  return *this;
}

Scatter &Scatter::Setabc()
{
  // See if the input name is in $CLIB/atomsf.lib and
  // set the form factors.
  
  string filename;
  
  if (getenv("CLIBD")  == NULL)
    Bp3Error("Scatter::Setabc","CCP4 environment variable $CLIBD not set");
  else
    filename            = getenv("CLIBD");

  filename             += "/atomsf.lib";
  
  std::ifstream atomsffile(filename.c_str(), std::ios::in);
  
  // If the input name is not in the above file, then exit.

  if (!atomsffile)
    Bp3Error("Scatter::Setabc",
	     "File " + filename + " could not be opened for reading");
  else
  {
    const int linesize  = 90;
    const int arraysize = 91;
    char line[arraysize];
    
    // Skip the first 32 lines of the library file to get to the first atom name    
    for (int i          = 0; i < 31; i++)
      atomsffile.getline(line, linesize);

    string atomname;
    
    do
    {   
      atomsffile >> atomname;
      
      for (unsigned i   = 0; i < atomname.size(); i++)
	atomname[i]     = toupper(atomname[i]);
      
      if (name         == atomname)
      {
	atomsffile.getline(line, linesize);
	
	double junk1, junk2;

	atomsffile >> junk1    >> junk2    >> c;
	atomsffile >> a[0] >> a[1] >> a[2] >> a[3];
	atomsffile >> b[0] >> b[1] >> b[2] >> b[3];
	atomsffile >> fp[0] >> fpp[0];
      }
      else
      {
	// We didn't find the atom we wanted, so go to the next atom name
	
	for (int i      = 0; i < 5; i++)
	  if (!atomsffile.getline(line, linesize))
	    Bp3Error("Scatter::Setabc",
		     "atom name " + name + " not found in $CLIBD/atomsf.lib");
      }
    } while (name      != atomname);
  }
  return *this;
}

void Scatter::print(const bool all) const
{
  printf("Scattering factors for atom %s\n", name.c_str());
  printf("A: %f %f %f %f\n", a[0], a[1], a[2], a[3]);
  printf("B: %f %f %f %f\n", b[0], b[1], b[2], b[3]);
  printf("C: %f\n", c);
  printf("Crystal #%u\n", xtal+1);
  for (unsigned w = 0; w < nwave; w++)
  {
    printf("Wavelength #%u\n", w+1);
    printf("fp = %7.3f fpp = %7.3f\n",fp[w], fpp[w]);
  }
  printf("\n");
  fflush(stdout);
}

// Site methods

Site::Site(const unsigned numin, const double xin, 
	   const double yin, const double zin, 
	   const bool refx, const bool refy, const bool refz)
{
  MINX = -THOUSAND; MAXX = THOUSAND;
  Setnumber(numin);
  Setx(0,xin).Setx(1,yin).Setx(2,zin);
  Setrefinex(0,refx).Setrefinex(1,refy).Setrefinex(2,refz);
}

vector<double> Site::ortho(vector<vector<work> > &frac2or) const
{
  // return site in orthogonal coordinates
  vector<double> orth(3, ZERO), fracx(3, ZERO);

  for (unsigned i   = 0; i < fracx.size(); i++)
    fracx[i]        = x[i] - (int)x[i];

  for (unsigned i   = 0; i < orth.size(); i++)
    for (unsigned j = 0; j < orth.size(); j++)
      orth[i]      += frac2or[i][j]*fracx[j];
  
  return orth;
}

void Site::invert(const unsigned sgnumber)
{
  if (sgnumber      == 81)   // I4(1)
  {
    x[0] = ONE      - x[0];
    x[1] = HALF     - x[1];
    x[2] = ONE      - x[2];    
  }
  else if (sgnumber == 98) // I4(1)22
  {
    x[0] = ONE      - x[0];
    x[1] = HALF     - x[1];
    x[2] = QUARTER  - x[2];
  }
  else if (sgnumber == 210) // F4(1)32
  {
    x[0] = QUARTER  - x[0];
    x[1] = QUARTER  - x[1];
    x[2] = QUARTER  - x[2];
  }
  else
    for (unsigned i = 0; i < 3; i++) x[i] = -x[i];
}

void Site::print() const
{
  printf("Site # %u\n", number+1);
  /*
  string refined( (refinex[0]) ? "X " : "" );
  refined.append( (refinex[1]) ? "Y " : "" );
  refined.append( (refinex[2]) ? "Z " : "" );
  printf(( (refined != "") ? "Refine " : "X Y Z Fixed ") << refined << std::endl;
  */
  printf("     X      Y      Z\n");
  printf("  %5.3f  %5.3f  %5.3f\n",x[0], x[1], x[2]);
  fflush(stdout);
}

// Atom methods

Atom::Atom()
{name  = "";   nsite   = 0; nform      = 0;    crystal   = 0;
 MINO  = ZERO; MAXO    = 50.0; MINB    = 2.500;  MAXB    = 250.0;
 MINAB = -250.0; MAXAB = 250.0; number = ONE; refinenumb = false;
}

Atom::Atom(const string namein, const unsigned nsitein, const unsigned nformin, 
	   const unsigned crysin, const bool isoin, const vector<double> &bfacin,
	   const double occin, const bool refbfacin, const bool refoccin)
{
 MINO                      = ZERO; MAXO    = 50.0; MINB    = 2.500; MAXB    = 250.0;
 MINAB                     = -250.0; MAXAB = 250.0; number = ONE;

  refinenumb                = false;
  Setname(namein);

  nsite                     = nsitein;
  nform                     = nformin;
  crystal                   = crysin;

  isotropic                 = isoin;
  
  if (isotropic)
    Setbiso(bfacin[0]);
  else
  {
    if (bfacin.size()      == 6)
      for (unsigned b       = 0; b < 6; b++)
	uaniso[b]           = bfacin[b];
    else if (bfacin.size() == 1)
      convertuaniso(bfacin[0]);
    else
      Bp3Error("Atom::Atom", "Incorrect vector size for anisotropic B-factor");
  }
   
  Setocc(occin).Setrefinebfac(refbfacin).Setrefineocc(refoccin);
}

bool Atom::posdefu()
{
  // check to see if the U matrix is positive (semi) definite
  
  return ((uaniso[0] >= ZERO) && ((uaniso[0]*uaniso[3] - uaniso[1]*uaniso[1]) >= ZERO) &&
	  ((-uaniso[2]*uaniso[2]*uaniso[3] + TWO*uaniso[1]*uaniso[2]*uaniso[4]
	    -uaniso[0]*uaniso[4]*uaniso[4] - uaniso[1]*uaniso[1]*uaniso[4]
	    +uaniso[0]*uaniso[3]*uaniso[5]) >= ZERO));
}

Atom &Atom::convertuaniso(const double value)
{
  // Set anisotropic B-factor from an isotropic one
  // ***NSP
  
  double val = value/(FOUR*TWOPI2);
  uaniso[0]  = between(MINAB, val, MAXAB);
  uaniso[1]  = ZERO;
  uaniso[2]  = ZERO;
  uaniso[3]  = between(MINAB, val, MAXAB);
  uaniso[4]  = ZERO;
  uaniso[5]  = between(MINAB, val, MAXAB);

  return *this;
}

void Atom::print() const
{
  /*
  string refined( (refineocc)  ? "O " : "" );
  refined.append( (refinebfac) ? "B " : "" );
  
  printf(( (refined != "") ? "Refine " : "O B Fixed ") << refined << std::endl;
  */

  if (number != ONE)
    printf("\nNumber of expected atoms: %6.2f\n",number);

  if (isotropic)
  {
    printf("     O        B\n");
    printf("  %7.3f  %7.3f\n",occ,biso);
  }
  else
  {
    printf("     O        U11      U22      U33      U12      U13      U23 \n");
    printf("  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f\n", occ, 
	   uaniso[0], uaniso[3], uaniso[5], uaniso[1], uaniso[2], uaniso[4]);
  }
  printf("\n");
  fflush(stdout);
}

// Model methods

void Model::inverthand(const unsigned sgnumber)
{
  printf("Inverting hand of coordinates\n\n");
  
  for (unsigned s   = 0; s < site.size(); s++)
    site[s].invert(sgnumber);

  print();
}

void Model::print(const bool all) const
{

  for (unsigned a   = 0; a < atom.size(); a++)
  {
    printf("Atom  %s\n", atom[a].Getname().c_str());
    if (site.size() > 0)
      site[atom[a].Getnsite()].print();
    atom[a].print();
  }

  for (unsigned f   = 0; f < form.size(); f++)
  {
    bool refine     = false;
    for (unsigned w = 0; w < form[f].Getnwave(); w++)
      if (form[f].Getrefinefp(w) || form[f].Getrefinefpp(w))
	refine      = true;
    if (refine     || all)
      form[f].print(all);
  }
}
