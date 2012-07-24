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

#include "sfdata.h"

Sfdata::Sfdata(const unsigned mbin, const unsigned refbin)
{
  MINK      = ZERO;  MAXK    = 100.0; MINBSCLISO = -200.0;  MAXBSCLISO = 200.0;
  MINDLUZ   = ZERO;  MAXDLUZ = 25.0;  MINADLUZ   = ZERO;    MAXADLUZ   = ONE;
  MINELUZ   = ZERO;  MAXELUZ = 25.0;  MINSDLUZ   = ZERO;    MAXSDLUZ   = 25.0;
  MINBINS   = 8;     MAXBINS = mbin;  MINRES     = ZERO;    MAXRES     = 1000.0;
  REFBINS   = refbin;
  defsdluz  = 0.65;
  anomalous = false; noniso  = true;
  nwave     = nxtal          = nbins             = nref                = 0;
  refined   = refinead       = refinesd          = true;
  refinee   = refinedmod     = false;
  kscale    = ONE;
  biso      = ZERO;
  // baniso    = {ZERO, ZERO, ZERO, ZERO, ZERO, ZERO};
  refinek   = refineb        = true;
  aniso     = false;
  hires     = MAXRES;
  lowres    = MINRES;
  sigmacut  = anosigmacut    = isosigmacut       = -1000000.1;
  anormscut = 3.00;
  sdatap    = sdevp          = sdatam            = sdevm               = "";
  type      = "";
}

void Sfdata::resize(const unsigned maxref, const bool anom)
{
  anomalous                    = anom;
  if (datap.size()            != maxref)
  {
    datap.resize(maxref, NOTUSED);
    devp.resize(maxref, NOTUSED);
    fcalcp.resize(maxref, ZERO);
    pcalcp.resize(maxref, ZERO);

    if (sfmodel.size())
    {
      fmodel.resize(maxref);
      pmodel.resize(maxref);
      // ***NSP
      refinedmod               = true;

      for (unsigned r          = 0; r < maxref; r++)
      {
	fmodel[r].resize(sfmodel.size(), ZERO);
	pmodel[r].resize(sfmodel.size(), ZERO);
      }
    }

    if (sfom.size())
      fom.resize(maxref, ZERO);
    if (sphi.size())
      phib.resize(maxref, ZERO);

    if (shla.size())
      hla.resize(maxref, ZERO);
    if (shlb.size())
      hlb.resize(maxref, ZERO);
    if (shlc.size())
      hlc.resize(maxref, ZERO);
    if (shld.size())
      hld.resize(maxref, ZERO);
  }
  
  if (anomalous)
  {
    if (datam.size()          != maxref)
    {
      datam.resize(maxref, NOTUSED);
      devm.resize(maxref, NOTUSED);
      fcalcm.resize(maxref, ZERO);
      pcalcm.resize(maxref, ZERO);
    }
  }
  else // // if there is no anomalous data, don't refine adluz!!
  {
    refinead                   = false;
    refinesd                   = false;
  }
}

void Sfdata::resizeclipper()
{
  // store only clipper objects - not bp3 ones..

  if (datap.size())
  {
    datap.resize(0);
    devp.resize(0);
    fcalcp.resize(0);
    pcalcp.resize(0);
    if (sfmodel.size())
    {
      fmodel.resize(0);
      pmodel.resize(0);
    }
  }
  
  if (anomalous)
    if (datam.size())
    {
      datam.resize(0);
      devm.resize(0);
      fcalcm.resize(0);
      pcalcm.resize(0);
    }
}

void Sfdata::resizebin()
{
  if (sigman.size() != nbins)
  {
    sigman.resize(nbins, ZERO);
    sigmanref.resize(nbins, ZERO);
    sigmadano.resize(nbins, ZERO);
    sigmap.resize(nbins, ZERO);
    sigmah.resize(nbins, ZERO);
    nshl.resize(nbins, 0);
    anonshl.resize(nbins, 0);
    astolsq.resize(nbins, ZERO);
    fovers.resize(nbins, ZERO);
    danoovers.resize(nbins, ZERO);  
  }
}

Sfdata &Sfdata::Setdluz(const unsigned s, const double din)
{
  if (s     < nbins)
    dluz[s] = between(MINDLUZ,din,MAXDLUZ);
  else
    Bp3Error("Sfdata::Setdluz", "out of bounds array indices");
  
  return *this;
}

Sfdata &Sfdata::Setdmod(const unsigned s, const unsigned m, const double din)
{
  if ( (s      < nbins) && (m < sfmodel.size()) )
    dmod[s][m] = between(MINDLUZ,din,MAXDLUZ);
  else
    Bp3Error("Sfdata::Setdmod", "out of bounds array indices");
  
  return *this;
}

Sfdata &Sfdata::Seteluz(const unsigned s, const double ein)
{
  if (s     < nbins )
    eluz[s] = between(MINELUZ,ein,MAXELUZ);
  else
    Bp3Error("Sfdata::Seteluz", "out of bounds array indices");
  
  return *this;
}

Sfdata &Sfdata::Setadluz(const unsigned s, const double adin)
{
  if (s      < nbins )
    adluz[s] = between(MINADLUZ,adin,MAXADLUZ);
  else
    Bp3Error("Crystal::Setadluz", "out of bounds array indices");
  
  return *this;
}

Sfdata &Sfdata::Setsdluz(const unsigned s, const double sdin)
{
  if (s      < nbins )
    sdluz[s] = between(MINSDLUZ,sdin,MAXSDLUZ);
  else
    Bp3Error("Crystal::Setsdluz", "out of bounds array indices");
  
  return *this;
}

Sfdata &Sfdata::Setdluz(const vector<double> &din)
{
  for (unsigned i = 0; i < din.size(); i++)
    dluz.push_back(between(MINDLUZ,din[i],MAXDLUZ));
  
  return *this;
}

Sfdata &Sfdata::Seteluz(const vector<double> &ein)
{
  for (unsigned i = 0; i < ein.size(); i++)
    eluz.push_back(between(MINELUZ,ein[i],MAXELUZ));
  
  return *this;
}

Sfdata &Sfdata::Setadluz(const vector<double> &adin)
{
  for (unsigned i = 0; i < adin.size(); i++)
    adluz.push_back(between(MINADLUZ,adin[i],MAXADLUZ));
  
  return *this;
}

Sfdata &Sfdata::Setsdluz(const vector<double> &sdin)
{
  for (unsigned i = 0; i < sdin.size(); i++)
    sdluz.push_back(between(MINSDLUZ,sdin[i],MAXSDLUZ));
  
  return *this;
}

/*
Sfdata &Sfdata::Setsdluzref(const vector<double> &sdin)
{
  for (unsigned i = 0; i < sdin.size(); i++)
    sdluzref.push_back(between(MINSDLUZ,sdin[i],MAXSDLUZ));
  
  return *this;
}
*/

Sfdata &Sfdata::Setnref()
{
  nref            = 0;
  for (unsigned r = 0; r < datap.size(); r++)    
    if (use(r))
      nref++;

  printf("Number of reflections in for Dataset %s is %u\n\n",name.c_str(), nref);
  
  return *this;
}

Sfdata &Sfdata::Setshell()
{
  if (lowres > EPSILON)
    lowsq    = QUARTER/(lowres*lowres);
  else
    Bp3Error("Sfdata::Setshell", "low resolution limit is zero");
  
  if (hires  > EPSILON)
    hisq     = QUARTER/(hires*hires);
  else
    Bp3Error("Sfdata::Setshell","high resolution limit is zero");

  if (hisq  != lowsq)
    vshell   = (work)(nbins)/(hisq - lowsq);
  else
    Bp3Error("Sfdata::Setshell","high and low resolution limits are the same");

  return *this;  
}

void Sfdata::checkparameters()
{
  // do not want to refine reference set scale and Luzzati parameters
  if ((!nxtal) && (!nwave))
  {
    refinek           = false;
    refineb           = false;
    refined           = false;
    refinesdref       = false;
  }

  /*     ***NSP
  if ( (dluz.size()  != nbins) || (adluz.size() != nbins) || 
       (eluz.size()  != nbins) || (sdluz.size() != nbins)   )    
    Bp3Error("Sfdata::checkparameters",
	     "nbins and Luzzati parameters don't match");
  */
}

void Sfdata::printLuzzati(const string target, const int cyc, const bool verbose) const
{
  if (verbose)
  {
    if (cyc           > 0)
      printf(" $TABLE: Cycle %d Dataset %s Luzzati parameters vs. STOL2 - iso, ano:\n",
	     cyc, name.c_str());
    else
      printf(" $TABLE: Dataset %s Luzzati parameters vs. STOL2 - iso, ano:\n",
	     name.c_str());
    
    printf("$GRAPHS");
  }
  
  if ( (target        == "UNCO") || (target == "") )
    printf(": Luzzati iso parameters vs. Res:N:4,5:\n");

  if (anomalous)
    if ( (target      == "UNCO") || (target == "") )
      printf(": Luzzati anom parameters vs. Res:N:4,6:\n");
    else
      printf(": Luzzati anom parameters vs. Res:N:4,5:\n");

  if ( target         == "MLHL")
    printf (": Luzzati model parameters vs. Res:N:4,5:\n");
      
  printf("$$\n Bin   LoRes  HiRes     Res ");
  if ( (target        == "UNCO") || (target == "") )
    printf("  Isomorphous");
  if (anomalous)
    printf("   Anomalous  $$\n$$\n");
  else
    printf("$$\n$$\n");
  for (unsigned s      = 0; s < nbins; s++)
  {
    printf(" %2d  %6.2f  %6.2f  %7.3f ",s+1,Getlores(s),Gethires(s),Getavres(s));
    if ( (target      == "UNCO") || (target == "") )
      printf("%10.5f  ",dluz[s]);
    else if (dmod.size())
      printf("%10.5f  ",dmod[s][0]);
    if (anomalous)
    {
      if ( (target   == "UNCO") || (target == "") )
	printf("%10.5f     ", adluz[s]);
      else if ( (target == "SAD") || (target == "MAD") )
	printf("%10.5f    ", sdluz[s]);
    }
    printf("\n");
  }
  printf("$$\n\n\n");
  fflush(stdout);  
}

void Sfdata::printnorm(const bool verbose) const
{
  printf("\n");
  if (verbose)
  {
    printf(" $TABLE: Dataset %s Normalization statistics vs. Resolution:\n", name.c_str());
    printf("$GRAPHS: SigmaN    vs. Res:N:4,6:\n");
  
    printf("       : SigmaH    vs. STOL2:N:4,8:\n");
    if (anomalous)
      printf("       : SigmaDano vs. STOL2:N:4,9:\n");
    printf("       : Variance  vs. STOL2:N:4,10:\n");
  }
  
  printf("$$\n Bin  LoRes  HiRes   Res    STOL2     Refls    SigmaN    SigmaNref    ");
  printf("  SigmaH     SigmaP     SigmaDano   Variance $$\n$$\n");
  for (unsigned s  = 0; s < nbins; s++)
  {
    if (anomalous)
    {
      printf(    " %2d  %6.2f %6.2f %6.2f  %7.5f   %5i  %10.2f  %10.2f  %10.2f  %10.2f  %10.2f",
		 s+1, Getlores(s), Gethires(s), Getavres(s), astolsq[s], nshl[s],
		 sigman[s], sigmanref[s], sigmah[s], sigmap[s], sigmadano[s]);
      printf(" %10.2f\n", ( (sigman[s] - dluz[s]*dluz[s]*sigmanref[s]) < EPSILON ) ?
	     0.0 : sigman[s] - dluz[s]*dluz[s]*sigmanref[s]);
    }
    else
    {
      printf(      " %2d  %6.2f %6.2f %6.2f  %7.5f   %5i  %10.2f  %10.2f  %10.2f  %10.2f  %10.2f",
		   s+1, Getlores(s), Gethires(s), Getavres(s), astolsq[s], nshl[s],
		   sigman[s], sigmanref[s], sigmah[s], sigmap[s], 0.00);
      printf(" %10.2f\n", ( sigman[s] - dluz[s]*dluz[s]*sigmanref[s]-sigmah[s] < EPSILON ) ?
	     ZERO : sigman[s] - dluz[s]*dluz[s]*sigmanref[s]-sigmah[s]);
    }
  }
  printf("$$\n\n");
  fflush(stdout);  
}

void Sfdata::printsignal() const
{
  printf("\n");
  printf(" $TABLE: Dataset %s Signal to noise ratio vs. Resolution:\n", name.c_str());
  string ston( ( (type == "AMPLITUDE") ? "F/SigmaF" : "I/SigmaI" ) );
  printf("$GRAPHS: %s    vs. Res:N:4,7:\n",ston.c_str());
  printf("       : %s    vs. STOL2:N:5,7:\n",ston.c_str());
  if (anomalous)
  {
    printf("       : Dano/SigDano vs. Res:N:4,9:\n");
    printf("       : Dano/SigDano vs. STOL2:N:5,9:\n");
    
  }
  

  printf("$$\n Bin  LoRes  HiRes    Res   STOL2    Refls     %s",ston.c_str());
  if (anomalous)
    printf("   AnoRefls    Dano/SigDano");
  printf(" $$\n$$\n");
  double total(ZERO), totalsignal(ZERO), totalano(ZERO), totalanosig(ZERO);
  for (unsigned s  = 0; s < nbins; s++)
  {
    printf(" %2d  %6.2f %6.2f %6.2f  %7.5f   %5i  %10.2f",
	   s+1, Getlores(s), Gethires(s), Getavres(s), astolsq[s], nshl[s],
	   fovers[s]);
    total         += (double)nshl[s];
    totalsignal   += ((double)nshl[s])*fovers[s];
    if (anomalous)
    {
      printf("     %5i    %10.2f\n",anonshl[s], danoovers[s]);
      totalano    += (double)anonshl[s];
      totalanosig += ((double)anonshl[s])*danoovers[s];
    }
    else
      printf("\n");
  }
  printf("$$\nTOTAL\n");
  printf("     %6.2f %6.2f                   %5.0f  %10.2f",
	 Getlores(0),Gethires(nbins-1), total, totalsignal/total);
	 //	 Getlores(0),Gethires(nbins-1), total, ( (total > ZERO) ? totalsignal/total : ZERO));
  if (anomalous)
    printf("     %5.0f    %10.2f", totalano, totalanosig/totalano);
    
    // printf("%5d  %7.4f", totalano, ( (totalano > ZERO) ? totalanosig/totalano : ZERO) );
  printf("\n\n");
  fflush(stdout);  
}

void Sfdata::print() const
{
  printf("Structure factor data\n");
  printf("Number of reflections: %u\n", nref);
  fflush(stdout);  
}
