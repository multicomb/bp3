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
#include "bp3likelihood.h"
#include "mytimer.h"
extern Timer Tsadgradient;

double Bp3likelihood::uncorrelatedgradient(const bool check, const bool outputmtz)
{
  // Calculates likelihood assuming uncorrelated errors.

  likelihood                               = ZERO;

  double REJECT((check && !outputmtz) ? DSMALL : SMALLESTD);
  
  const bool deriv(true);
  
  unsigned flagged(0), nstop(0);
  double stopfom(ZERO);
  unsigned dmax(xtal.sf.size());

  // use only datasets from one crystal in phasing for mad
  if (xtal.mad && xtal.dataset.size()      > 1)
      dmax                                 = xtal.dataset[1][0];

  if (check)
    xtal.Setselected( (outputmtz ? "PHASE" : "REFINE") );

  CMtz::MTZ *MTZIN(NULL), *MTZOUT(NULL);

  if (outputmtz)
  {
    MTZOUT                                 = CMtz::MtzMalloc(0,0);
    if (allin)
    {
      if (getenv("HKLIN")                 != NULL)
	MTZIN                              = CMtz::MtzGet("HKLIN",0);
      else if (xtal.sf[0].mtzin.size())
	MTZIN                              = CMtz::MtzGet(xtal.sf[0].mtzin.c_str(),0);

      if (!MTZIN)
	Bp3Error("Bp3likelihood::uncorrelatedgradient","Can not read input mtz");
    }
    
    setupmtz(MTZOUT,MTZIN);
  }

  unsigned inc((allin) ? colin.size() : 3);
  
  for (unsigned d                          = 0; d < dmax  ; d++)
  {
    initstat();
    for (unsigned s                        = 0; s < xtal.sf[d].nbins; s++)
    {
      dLdsigmah[d][s]                      = dLdsumfpp[d][s] = ZERO;
      dLddluz[d][s]                        = dLdadluz[d][s]  = ZERO;
      anshl[d][s]                          = cnshl[d][s]     = 0;
      if (d                               == 0)
	afom[s]                            = cfom[s]         = ZERO;
      sanonshl[d][s]                       = 0;
    }
  }

  vector<double> koverall(dmax, ONE), var(dmax, ZERO);
  vector<double> besselarg(dmax), avar(dmax);
  vector<double> fcalc(dmax), fcalcp(dmax), fcalcm(dmax);
  vector<double> diff(dmax, ZERO), adiff(dmax);
  vector<double> tsim(dmax, ZERO);
  vector<double> cosp(dmax), sinp(dmax);
  vector<double> cosm(dmax), sinm(dmax);
  vector<double> cosharg(dmax);
  vector<double> fcalc2(dmax), diff2(dmax, ZERO);
  vector<double> nfcalc1(dmax), nfcalc2(dmax);
  
  vector<double> dLdd(dmax, ZERO),  dLdvar(dmax, ZERO);
  vector<double> dLdad(dmax, ZERO), dLdavar(dmax, ZERO);
  vector<double> dLdFp(dmax, ZERO), dLdPp(dmax, ZERO);
  vector<double> dLdFm(dmax, ZERO), dLdPm(dmax, ZERO);

  vector<double> dphasedfp(dmax, ZERO), dphasedpp(dmax, ZERO);
  vector<double> dphasedfm(dmax, ZERO), dphasedpm(dmax, ZERO);
  vector<double> dphaseddluz(dmax, ZERO), dphasedvar(dmax, ZERO);
  vector<double> dphasedadluz(dmax, ZERO), dphasedavar(dmax, ZERO);
  vector<double> delcosp(dmax), delcosm(dmax);

  vector<double> isolofc(dmax, ZERO), alofc(dmax, ZERO);
  vector<double> iphaselofc(dmax, ZERO), aphaselofc(dmax, ZERO);
  vector<double> ecosdiff(dmax, ZERO), ephasecosdiff(dmax, ZERO);
  
  for (unsigned r                          = 0; r < xtal.maxselref; r++)
  {
    // Default value for columns to be written out is MNF
    vector<float> fdata(9+inc, CCP4::ccp4_nan().f);

    if (outputmtz)
      storeinitialcolumns(MTZIN, fdata, inc, r);

    double eps((double) xtal.epsilon[r]);
    
    unsigned counts(0);
    int mders(-1);
    for (unsigned d                        = 0; d < dmax; d++)
      if (xtal.sf[d].use(r))
      {
 	counts++;
	mders++;
	if (xtal.sf[d].anomalous          && !xtal.centric[r])
	  if (xtal.sf[d].anouse(r))
	    counts++;
      }
    
    if (outputmtz                         && (counts == 1))
      for (unsigned d                      = 0; d < dmax; d++)
	if (xtal.sf[d].use(r))
	  counts                           = (xtal.sf[d].sigmah[0] > ZERO) ? 1 : 0;
    
    if ( (counts                           > 1) && xtal.selected[r])
    {
      unsigned sd(0);
      for (unsigned d                      = 0; ( ((sd=d) < dmax) && 
						  !(xtal.sf[d].use(r)) ); d++); 

      double fos(xtal.sf[sd].foversigma(r));

      bool onedim(xtal.sf[0].use(r) && (fos > threshold) && !xtal.heavyref);

      unsigned ind((onedim) ? 1 : sd);
      
      for (unsigned d                      = 0; d < dmax; d++)
      {
	dLdAp[d][r]                        = dLdBp[d][r] = ZERO;
	if (xtal.sf[d].anomalous)
	  dLdAm[d][r]                      = dLdBm[d][r] = ZERO;

	dLdFp[d]                           = dLdPp[d]  = dLdFm[d] = dLdPm[d]   = ZERO;
	dLdd[d]                            = dLdvar[d] = dLdad[d] = dLdavar[d] = ZERO;

	if (xtal.sf[d].use(r)  || (d      == 0) )
	  if (xtal.centric[r])
	    cnshl[d][xtal.bin(d,r)]++;
	  else
	    anshl[d][xtal.bin(d,r)]++;

	isolofc[d]                         = alofc[d]  = ecosdiff[d]           = ZERO;
      }
      double totcos(ZERO), totsin(ZERO), totamp(ZERO), totampcos(ZERO), totampsin(ZERO);
      double totlogcos(ZERO), totlogsin(ZERO), totlogcos2(ZERO), totlogsin2(ZERO);
      
      double integral(ZERO), product(ONE);
      vector<double> dervar(dmax, ZERO);
      vector<double> weight, node;

      if (!(xtal.centric[r] || xtal.onlycentrics))
      {
	if (onedim)
	{
	  node.push_back(ZERO);
	  weight.push_back(ONE);
	}
	else
	{
	  weight                           = afweight;
	  node                             = afnode;
	}
	
	vector<bool> deravar(dmax, true);

	for (unsigned d                    = ind; d < dmax; d++)
	  if (xtal.sf[d].use(r))
	  {
	    unsigned s(xtal.bin(d,r));
	    var[d]                         = pow(koverall[d]*xtal.sf[d].devmean(r), 2);
	    
	    if (xtal.sf[d].noniso         && (d != sd) )
	    {
	      double temp(eps*xtal.sf[d].var(koverall[d],s,xtal.heavyref));
	      dervar[d]                    = ONE;
	      if (temp                     > ZERO)
                var[d]                    += temp;
	      else
		dervar[d]                  = ZERO;
	    }
 	    besselarg[d]                   = TWO*koverall[d]*xtal.sf[d].datamean(r)/var[d];
	    cosp[d]                        = tab.Cos(xtal.sf[d].pcalcp[r]);
	    sinp[d]                        = tab.Sin_charged(xtal.sf[d].pcalcp[r]);
	    if (mders)
	      product                     *= besselarg[d];
	    asumfph[d][s]                 += koverall[d]*xtal.sf[d].datamean(r);
	    asumfh[d][s]                  += xtal.sf[d].fcalcp[r];
	    if (xtal.sf[0].use(r))
	      asumdiff[d][s]              += fabs(koverall[d]*xtal.sf[d].datamean(r) -
						  koverall[0]*xtal.sf[0].datamean(r));
 	    if (xtal.sf[d].anomalous)
	    {
	      cosm[d]                      = tab.Cos(xtal.sf[d].pcalcm[r]);
	      sinm[d]                      = tab.Sin_charged(xtal.sf[d].pcalcm[r]);
	      if (xtal.sf[d].anouse(r))
	      {
		avar[d]                    = std::max(eps*(koverall[d]*koverall[d]*
							   xtal.sf[d].sigmadano[s]*
							   (ONE-xtal.sf[d].adluz[s]*
							    xtal.sf[d].adluz[s])),ZERO);
		deravar[d]                 = (avar[d] > ZERO);
		avar[d]                   += pow(koverall[d]*xtal.sf[d].sigdano(r), 2);
		product                   /= sqrt(PI*avar[d]);
		sumimfh[d][s]             += fabs(xtal.sf[d].fcalcp[r]*
						  tab.Sin_charged(xtal.sf[d].pcalcp[r]));
		sumdano[d][s]             += fabs(xtal.sf[d].dano(r));
		anosumfph[d][s]           += fabs(xtal.sf[d].datamean(r));
		sanonshl[d][s]++;
	      }
	    }
	  }

	for (unsigned i                    = 0; i < pweight.size(); i++)
	{
	  double afvarnode(ZERO), afvarnode1(ZERO), saf1(ONE);
	  double phasesum(ZERO), phaseamp(ZERO);

	  for (unsigned d                  = ind; d < dmax; d++)
 	    if (xtal.sf[d].use(r))
	    {
	      delcosp[d]                   = gcos[i]*cosp[d] + gsin[i]*sinp[d];
	      if (xtal.sf[d].anomalous)
		delcosm[d]                 = gcos[i]*cosm[d] + gsin[i]*sinm[d];
	      dphasedfp[d]                 = dphasedpp[d]   = dphasedfm[d]     = ZERO;
	      dphasedpm[d]                 = dphaseddluz[d] = dphasedvar[d]    = ZERO;
	      dphasedadluz[d]              = dphasedavar[d] = ZERO;
	      
	      iphaselofc[d]                = aphaselofc[d]  = ephasecosdiff[d] = ZERO;
	    }

	  for (unsigned a                  = 0; a < weight.size(); a++)
	  { 
 	    double expsum(ZERO), prob(ONE);

	    afvarnode                      = (xtal.sf[sd].devmean(r)*node[a] +
					      xtal.sf[sd].datamean(r));
	    if (xtal.heavyref)
	    {
	      afvarnode1                   = (afvarnode*afvarnode -
					      (xtal.sf[sd].fcalcp[r]*
					       xtal.sf[sd].fcalcp[r]*
					       (ONE - delcosp[sd]*
						delcosp[sd])));
	      saf1                         = sign(afvarnode1);
	      afvarnode1                   = std::max(sqrt(fabs(afvarnode1)), DSMALL);
	      afvarnode                    = (afvarnode1 -
					      xtal.sf[sd].fcalcp[r]*delcosp[sd]);
 	    }

	    for (unsigned d                = ind; d < dmax; d++)
	      if (xtal.sf[d].use(r))
	      {
		double df(xtal.sf[d].dluz[xtal.bin(d,r)]*afvarnode);
		fcalc[d] = fcalcp[d]       = sqrt(std::max(df*df+
					 		   xtal.sf[d].fcalcp[r]*
							   xtal.sf[d].fcalcp[r]+
							   TWO*df*delcosp[d]*
							   xtal.sf[d].fcalcp[r],
							   EPSILON));
		if (xtal.sf[d].anomalous)
		{
		  fcalcm[d]                = sqrt(std::max(df*df+
				 			   xtal.sf[d].fcalcm[r]*
							   xtal.sf[d].fcalcm[r]+
							   TWO*df*delcosm[d]*
							   xtal.sf[d].fcalcm[r],
							   EPSILON));
		  if (xtal.sf[d].anouse(r))
		  {
		    fcalc[d]               = (fcalcp[d] + fcalcm[d])*HALF;
		    adiff[d]               = (koverall[d]*xtal.sf[d].dano(r) - 
 					      (fcalcp[d] - fcalcm[d]));
		    expsum                -= adiff[d]*adiff[d]/avar[d];
		  }
		  else if (xtal.sf[d].usem(r))
		    fcalc[d]               = fcalcm[d];  
		}
		if (mders)
		{
		  diff[d]                  = koverall[d]*xtal.sf[d].datamean(r) - fcalc[d];
		  expsum                  -= diff[d]*diff[d]/var[d];
		  double io(tab.I0e(besselarg[d]*fcalc[d]));
		  tsim[d]                  = tab.Sim_preI0(besselarg[d]*fcalc[d], io);
		  prob                    *= io;
		}
	      }
	    prob                          *= tab.ExpM(-expsum);

	    if (outputmtz                 && xtal.heavyref)
	      prob                        *= (fabs(xtal.sf[sd].devmean(r)*node[a] +
			 			   xtal.sf[sd].datamean(r))/afvarnode1);
	    double ppweight(prob*weight[a]);
	    phasesum                      += ppweight;
	    phaseamp                      += afvarnode*ppweight;

	    for (unsigned d                = ind; d < dmax; d++)
 	      if (xtal.sf[d].use(r))
	      {
		unsigned s(xtal.bin(d,r));
		double dfjcalc(ZERO), dadiffdfp(ZERO);
		double temp                = ((tsim[d]*koverall[d]*xtal.sf[d].datamean(r)
					       - fcalc[d])*ppweight/var[d]);

		if (xtal.sf[d].anomalous)
		{
		  double dfjcalcm(ZERO), dadiffdfm(ZERO);

		  if (xtal.sf[d].anouse(r))
		  {
		    if (mders)
		      dfjcalc = dfjcalcm   = temp/TWO;

		    dadiffdfp              = adiff[d]*ppweight/avar[d];
		    dadiffdfm              = adiff[d]*ppweight/avar[d];
		    dphasedavar[d]        += adiff[d]*adiff[d]*ppweight;
		  }
		  else if (xtal.sf[d].usem(r)  && mders)
		    dfjcalcm               = temp;
		  else if (mders)
		    dfjcalc                = temp;

		  dphaseddluz[d]          += ((xtal.sf[d].dluz[s]*afvarnode +
					       xtal.sf[d].fcalcm[r]*
					       delcosm[d])/fcalcm[d]*
					      (dfjcalcm - dadiffdfm))*afvarnode;
		  dphasedfm[d]            += ((xtal.sf[d].fcalcm[r]+xtal.sf[d].dluz[s]*
					       afvarnode*delcosm[d])/fcalcm[d]*
					      (dfjcalcm - dadiffdfm));
		  dphasedpm[d]            += ((dfjcalcm - dadiffdfm)/fcalcm[d]*
					      afvarnode);

		  if (xtal.heavyref)
		  {
		    double dfcdf((xtal.sf[d].dluz[s]*
				  (xtal.sf[d].dluz[s]*afvarnode +
				   xtal.sf[d].fcalcm[r]*delcosm[d])
				  )*(dfjcalcm - dadiffdfm)/fcalcm[d]);
		    if (afvarnode1         > DSMALL)
		    {
		      dphasedfp[sd]       -= ((xtal.sf[sd].fcalcp[r]*saf1/afvarnode1*
                                               (ONE - delcosp[sd]*delcosp[sd])
                                               + delcosp[sd])*dfcdf);
		      dphasedpp[sd]       += ((delcosp[sd]*xtal.sf[sd].fcalcp[r]*saf1
                                               /afvarnode1 - ONE)*dfcdf);
		    }
		    else
		    {
		      dphasedfp[sd]       -= delcosp[sd]*dfcdf;
		      dphasedpp[sd]       -= dfcdf;
		    }
 		  }
		}
		else if (mders)
		  dfjcalc                  = temp;

		dphasedvar[d]             += (diff[d]*diff[d]/var[d] +
					      besselarg[d]*fcalc[d]*
					      (ONE - tsim[d]))*ppweight;
		dphaseddluz[d]            += ((xtal.sf[d].dluz[s]*afvarnode +
					       xtal.sf[d].fcalcp[r]*delcosp[d])
					      /fcalcp[d]*
					      (dfjcalc + dadiffdfp))*afvarnode;
		dphasedfp[d]              += ((xtal.sf[d].fcalcp[r]+xtal.sf[d].dluz[s]*
					       afvarnode*delcosp[d])/fcalcp[d]*
					      (dfjcalc + dadiffdfp));
		dphasedpp[d]              += ((dfjcalc + dadiffdfp)/fcalcp[d]*
					      afvarnode);

		if (xtal.heavyref)
		{
		  double dfcdf(xtal.sf[d].dluz[s]*
			       (xtal.sf[d].dluz[s]*afvarnode +
				xtal.sf[d].fcalcp[r]*delcosp[d]
				)*(dfjcalc + dadiffdfp)/fcalcp[d]);
		  if (afvarnode1           > DSMALL)
		  {
		    dphasedfp[sd]         -= (xtal.sf[sd].fcalcp[r]*saf1/afvarnode1*
                                              (ONE - delcosp[sd]*delcosp[sd]) +
                                              delcosp[sd])*dfcdf;
		    dphasedpp[sd]         += (delcosp[sd]*xtal.sf[sd].fcalcp[r]*saf1
                                              /afvarnode1 - ONE)*dfcdf;
		  }
		  else
		  {
		    dphasedfp[sd]         -= delcosp[sd]*dfcdf;
		    dphasedpp[sd]         -= dfcdf;		    
		  }
		}
	      }

	    for (unsigned d                = ind; d < dmax; d++)
	      if (xtal.sf[d].use(r))
	      {
		ephasecosdiff[d]          += delcosp[d]*ppweight;
		iphaselofc[d]             += fabs(diff[d])*ppweight;
		if (xtal.sf[d].anomalous)
		  if (xtal.sf[d].anouse(r))
		    aphaselofc[d]         += fabs(adiff[d])*ppweight;  
	      }
	  }
	  integral                        += phasesum*pweight[i];
	  totcos                          += gcos[i]*phasesum*pweight[i];
	  totsin                          += gsin[i]*phasesum*pweight[i];
	  totamp                          += phaseamp*pweight[i];
	  totampcos                       += gcos[i]*phaseamp*pweight[i];
	  totampsin                       += gsin[i]*phaseamp*pweight[i];
	  if (phasesum                     > SMALLESTD)
	  {
	    double temp(log(phasesum)*pweight[i]);
	    totlogcos                     += gcos[i]*temp;
	    totlogsin                     += gsin[i]*temp;
	    totlogcos2                    += ((gcos[i]*gcos[i] - gsin[i]*gsin[i])*temp);
	    totlogsin2                    += (TWO*gcos[i]*gsin[i]*temp);
	  } 

	  for (unsigned d                  = ind; d < dmax; d++)
	    if (xtal.sf[d].use(r))
	    {
	      dLdd[d]                     += dphaseddluz[d]*pweight[i];
	      dLdvar[d]                   += dphasedvar[d]*pweight[i];
	      dLdFp[d]                    += dphasedfp[d]*pweight[i];
	      dLdPp[d]                    += (dphasedpp[d]*pweight[i]*
					      (gsin[i]*cosp[d] - gcos[i]*sinp[d]));
	      if (xtal.sf[d].anomalous)
	      {
		dLdavar[d]                += dphasedavar[d]*pweight[i];
		dLdFm[d]                  += dphasedfm[d]*pweight[i];
		dLdPm[d]                  += (dphasedpm[d]*pweight[i]*
					      (gsin[i]*cosm[d] - gcos[i]*sinm[d]));
	      }
	    }

	  for (unsigned d                  = ind; d < dmax; d++)
	    if (xtal.sf[d].use(r))
	    {
	      ecosdiff[d]                 += ephasecosdiff[d]*pweight[i];
	      isolofc[d]                  += iphaselofc[d]*pweight[i];
	      if (xtal.sf[d].anomalous)
		if (xtal.sf[d].anouse(r))
		  alofc[d]                += aphaselofc[d]*pweight[i];
	    }
	}
	
	if (integral                       > REJECT)
	  for (unsigned d                  = ind; d < dmax; d++)
	    if (xtal.sf[d].use(r))
	    {
	      unsigned s                   = xtal.bin(d,r);
	      double denom(ONE/(integral*var[d]));
	      dLdvar[d]                    = dervar[d]*(ONE/var[d] - dLdvar[d]*denom);
	      dLddluz[d][s]               -= TWO*dLdd[d]/integral;

	      if (mders                   && xtal.sf[d].noniso)
	      {
		dLddluz[d][s]             -= (TWO*eps*xtal.sf[d].dluz[s]*
					      dLdvar[d]*xtal.sf[d].sigmanref[s]);
		if (updatesigmah          && !xtal.heavyref)
		  dLdsigmah[d][s]         -= eps*dLdvar[d];
 	      }
	      denom                        = TWO/integral;
	      if (xtal.sf[d].anomalous)
	      {
		if (xtal.sf[d].anouse(r))
		{
		  dLdavar[d]               = (HALF - dLdavar[d]/(integral*avar[d])
					      )/avar[d];
		  if (deravar[d])
		    dLdadluz[d][s]        -= (xtal.sf[d].adluz[s]*dLdavar[d]*
					      eps*koverall[d]*koverall[d]*
					      xtal.sf[d].sigmadano[s])*TWO;
		}
		dLdAm[d][r]                = -(dLdFm[d]*cosm[d] -
				 	       dLdPm[d]*sinm[d]*xtal.sf[d].dluz[s])*denom;
		dLdBm[d][r]                = -(dLdFm[d]*sinm[d] +
					       dLdPm[d]*cosm[d]*xtal.sf[d].dluz[s])*denom;
	      }
	      dLdAp[d][r]                  = -(dLdFp[d]*cosp[d] -
					       dLdPp[d]*sinp[d]*xtal.sf[d].dluz[s])*denom;
	      dLdBp[d][r]                  = -(dLdFp[d]*sinp[d] +
					       dLdPp[d]*cosp[d]*xtal.sf[d].dluz[s])*denom;
	    }
	product                           *= integral*PI;

	if (!onedim)
	  product                         *= xtal.sf[sd].devmean(r);
      }
      else if (xtal.centric[r]            && !xtal.onlyacentrics)   // Centric case
      {
	if (onedim)
	{
	  node.push_back(ZERO);
	  weight.push_back(ONE);
	}
	else
	{
	  weight                           = cweight;
	  node                             = cnode;
	}

	vector<double> samephase(dmax, ONE);

	for (unsigned d                    = ind; d < dmax; d++)
	  if(xtal.sf[d].use(r))
	  {
	    unsigned s                     = xtal.bin(d,r);
	    var[d]                         = pow(koverall[d]*xtal.sf[d].devmean(r), 2);
	    if (xtal.sf[d].noniso         && (d != sd))
	    {
	      double temp(eps*xtal.sf[d].var(koverall[d],s,xtal.heavyref));
	      dervar[d]                    = (temp > ZERO);
	      if (temp                     > ZERO)
 		var[d]                    += temp;
	      else
		likelihood                += 100.0;
	    }
	    if (tab.Cos(xtal.centricphase[r] - xtal.sf[d].pcalcp[r]) < ZERO)
	      samephase[d]                 = -ONE;
	    product                       /= (TWO*PI*var[d]);
	    cosharg[d]                     = koverall[d]*xtal.sf[d].datamean(r)/var[d];
	    csumfph[d][s]                 += koverall[d]*xtal.sf[d].datamean(r);
	    csumfh[d][s]                  += xtal.sf[d].fcalcp[r];
	    if (xtal.sf[0].use(r))
	      csumdiff[d][s]              += fabs(koverall[d]*xtal.sf[d].datamean(r) -
	                                          koverall[0]*xtal.sf[0].datamean(r));
	  }

	for (unsigned i                    = 0; i < weight.size(); i++)
	{
	  double cvarnode(xtal.sf[sd].devmean(r)*node[i] +
		 	  xtal.sf[sd].datamean(r));
	  
	  double cvarnode1(cvarnode), cvarnode2(cvarnode);

	  if (xtal.heavyref)
	  {
	    cvarnode1                      = (cvarnode - samephase[sd]*xtal.sf[sd].fcalcp[r]);
	    cvarnode2                      = (cvarnode + samephase[sd]*xtal.sf[sd].fcalcp[r]);
	  }

	  double prob1(ONE),  expsum1(ZERO),   prob2(ONE),  expsum2(ZERO);
	  
	  for (unsigned d                  = ind; d < dmax; d++)
	    if (xtal.sf[d].use(r))
	    {
	      unsigned s(xtal.bin(d,r));	      
	      
	      nfcalc1[d]                   = (xtal.sf[d].dluz[s]*cvarnode1 +
					      samephase[d]*xtal.sf[d].fcalcp[r]);
	      nfcalc2[d]                   = (xtal.sf[d].dluz[s]*cvarnode2 -
					      samephase[d]*xtal.sf[d].fcalcp[r]);

	      fcalc[d]                     = fabs(nfcalc1[d]);
	      fcalc2[d]                    = fabs(nfcalc2[d]);

	      diff[d]                      = (koverall[d]*xtal.sf[d].datamean(r)
					      - fcalc[d]);
	      diff2[d]                     = (koverall[d]*xtal.sf[d].datamean(r)
					      - fcalc2[d]);
	      if (d                       != sd)
	      {
		expsum1                   -= diff[d]*diff[d]/(TWO*var[d]);
		expsum2                   -= diff2[d]*diff2[d]/(TWO*var[d]);
	      }
	      prob1                       *= ONE + tab.ExpM(TWO*cosharg[d]*fcalc[d]);
	      prob2                       *= ONE + tab.ExpM(TWO*cosharg[d]*fcalc2[d]);
	    }

	  prob1                           *= tab.ExpM(-expsum1);
	  prob2                           *= tab.ExpM(-expsum2);
	  
	  integral                        += (prob1 + prob2)*weight[i];
	  totcos                          += (prob1 - prob2)*weight[i];
	  totampcos                       += (prob1*cvarnode1 -
					      prob2*cvarnode2)*weight[i];
	  
	  for (unsigned d                  = ind; d < dmax; d++)
	    if (xtal.sf[d].use(r))
	      isolofc[d]                  += (prob1*fabs(diff[d]) + 
					      prob2*fabs(diff2[d]))*weight[i];

	  if (deriv)
	    for (unsigned d                = ind; d < dmax; d++)
	      if (xtal.sf[d].use(r))
	      {
		double tanharg1(tanh(cosharg[d]*fcalc[d]));
		double tanharg2(tanh(cosharg[d]*fcalc2[d]));

		double derivfjcalc1(ZERO), derivfjcalc2(ZERO);

		if (d                     != sd)
		{
		  derivfjcalc1             = prob1*(koverall[d]*xtal.sf[d].datamean(r)*
						    tanharg1 - fcalc[d])*weight[i];

		  derivfjcalc2             = prob2*(koverall[d]*xtal.sf[d].datamean(r)*
		 		 		    tanharg2 - fcalc2[d])*weight[i];
		}
		
		dLdvar[d]                 += (prob1*(diff[d]*diff[d]/var[d] +
						     TWO*fcalc[d]*koverall[d]*
						     xtal.sf[d].datamean(r)/var[d]*
						     (ONE - tanharg1)) +
					      prob2*(diff2[d]*diff2[d]/var[d] +
						     TWO*fcalc2[d]*koverall[d]*
						     xtal.sf[d].datamean(r)/var[d]*
						     (ONE - tanharg2)))*weight[i];
 		dLdd[d]                   += (sign(nfcalc1[d])*cvarnode1*derivfjcalc1 +
					      sign(nfcalc2[d])*cvarnode2*derivfjcalc2);

		double temp(sign(nfcalc1[d])*derivfjcalc1 -
			    sign(nfcalc2[d])*derivfjcalc2);
		  
  		dLdFp[d]                  += temp*samephase[d];
		
		if (xtal.heavyref)
		  dLdFp[sd]               -= xtal.sf[d].dluz[xtal.bin(d,r)]*temp*samephase[sd];
	      }
	  
	}

	product                            = sqrt(product)*integral;
	if (!onedim)
	  product                         *= xtal.sf[sd].devmean(r);
	
	if (deriv && (integral             > REJECT))
	  for (unsigned d                  = ind; d < dmax; d++)
	    if (xtal.sf[d].use(r))
	    {
	      unsigned s(xtal.bin(d,r));
	      double temp(ONE/(integral*var[d]));
	      dLdvar[d]                    = (ONE - dLdvar[d]/integral)*HALF;

	      dLddluz[d][s]               -= dLdd[d]*temp;

	      if (xtal.sf[d].noniso       && (d != sd) )
	      {
		  dLddluz[d][s]           -= (TWO*eps*xtal.sf[d].dluz[s]*dLdvar[d]*
					      xtal.sf[d].sigmanref[s]/var[d]);
		if (updatesigmah          && !xtal.heavyref)
		  dLdsigmah[d][s]         -= dLdvar[d]*eps/var[d];
	      }
	      dLdAp[d][r]                  = (-dLdFp[d]*tab.Cos(xtal.sf[d].pcalcp[r])
					      *temp);
	      dLdBp[d][r]                  = (-dLdFp[d]*tab.Sin_charged(xtal.sf[d].pcalcp[r])
					      *temp);
	    }
	}
      
      if (product                          > REJECT)
	likelihood                        -= log(product);
      else // A very low likelihood, so as not to confuse the minimizer
      {
	likelihood                        += 700.0;
	if (check)
	{
	  flagged++;
	  xtal.selected[r]                 = false;
	}
      }
    
      if (product                          > REJECT)
      {
	unsigned s0                        = xtal.bin(0,r);
	double fom                         = sqrt(totcos*totcos +
 						  totsin*totsin)/integral;
 	double fb                          = sqrt(totampcos*totampcos +
						  totampsin*totampsin)/integral;
	
	double phib(ZERO), hla(ZERO), hlb(ZERO);
	if (xtal.centric[r])
	{
	  phib                             = ((totcos >= ZERO) ? xtal.centricphase[r] :
					      xtal.centricphase[r] + PI);
	  cfom[s0]                        += fom;
	  for (unsigned d                  = 0; d < dmax  ; d++)
	  {
	    unsigned s                     = xtal.bin(d,r);
	    cisolofc[d][s]                += isolofc[d]/integral;
	  }
	  if (outputmtz)
	  {
	    double temp                    = atanh(fom);
	    hla                            = temp*tab.Cos(phib);
	    hlb                            = temp*tab.Sin_charged(phib);
	  }
	}
	else
	{
	  phib                             = atan2(totsin, totcos);
	  afom[s0]                        += fom;
	  for (unsigned d                  = 0; d < dmax  ; d++)
	  {
	    unsigned s                     = xtal.bin(d,r);
	    aisolofc[d][s]                += isolofc[d]/integral;  
	    if (xtal.sf[d].anomalous)
	      if (xtal.sf[d].anouse(r))
		anolofc[d][s]             += alofc[d]/integral;
	  }
	  if (outputmtz)
	  {
	    hla                            = totlogcos;
	    hlb                            = totlogsin;
	  }
	}
	stopfom                           += fom;
	nstop++;
	
	if (outputmtz)
	{
	  fdata[inc + 2]                   = (float) fb;
	  fdata[inc + 3]                   = (float) (phib/DEGREEtoRAD);
	  fdata[inc + 4]                   = (float) fom;
	  fdata[inc + 5]                   = (float) hla;
	  fdata[inc + 6]                   = (float) hlb;   
	  fdata[inc + 7]                   = (float) totlogcos2;
	  fdata[inc + 8]                   = (float) totlogsin2;
	}
      }
    }
    else if (outputmtz && (counts         == 1) && xtal.selected[r])
    {
      for (unsigned d                      = 0; d < dmax; d++)
	if (xtal.sf[d].use(r))
	{
	  double fobs, arg(ZERO), fom(CCP4::ccp4_nan().f);
	  unsigned s                       = xtal.bin(d,r);
	  double var                       = std::max(eps*(xtal.sf[d].sigman[s] - 
		   					   xtal.sf[d].sigmah[s]), EPSILON);
   	  if (xtal.sf[d].usep(r))
	  {
	    fobs                           = xtal.sf[d].datap[r];
	    var                           += xtal.sf[d].devp[r]*xtal.sf[d].devp[r];
	    arg                            = TWO*fobs*xtal.sf[d].fcalcp[r]/var;
	  }
	  else
	  {
	    fobs                           = xtal.sf[d].datam[r];
	    var                           += xtal.sf[d].devm[r]*xtal.sf[d].devm[r];
	    arg                            = TWO*fobs*xtal.sf[d].fcalcm[r]/var;
	  }

	  if (xtal.centric[r])
	  {
	    fom                            = tanh(arg/TWO);
	    cfom[s]                       += fom;
	    cnshl[d][s]++;
	  } 	
	  else
	  {
	    fom                            = sim(arg,0,0);
	    afom[s]                       += fom;
	    anshl[d][s]++;
	  }
	  
	  double phib                      = xtal.sf[d].pcalcp[r];

	  fdata[inc + 2]                   = (float) fom*fobs;
	  fdata[inc + 3]                   = (float) (phib/DEGREEtoRAD);
	  fdata[inc + 4]                   = (float) fom;
	  double temp                      = atanh(fom);
	  fdata[inc + 5]                   = (float) temp*cos(phib);
	  fdata[inc + 6]                   = (float) temp*sin(phib);
	  fdata[inc + 7]                   = (float) ZERO;
	  fdata[inc + 8]                   = (float) ZERO;
	}
    }
    if (outputmtz)
      CMtz::ccp4_lwrefl(MTZOUT, &fdata[0], &colout[0], fdata.size(), r+1);
  }

  if (outputmtz)
  {
    if (allin)
      CMtz::MtzFree(MTZIN);

    if (nstop)
    {
      char result[125];
      sprintf(result,"The overall FOM is %.3f\n\n",stopfom/((double)nstop));
      Bp3Result(result);
    }
    
    CMtz::MtzPut(MTZOUT, " ");
    CMtz::MtzFree(MTZOUT);
  }  
  
  if (xtal.verbose > 1)
    printf("LIKELIHOOD = %f\n", likelihood);
  
  if (check && flagged)
    printf("Number of reflections with low likelihood values: %u\n", flagged);
  
  return likelihood;
}

double Bp3likelihood::sadgradient(const bool check, const bool outputmtz)
{
  Tsadgradient.start("Main");
  // Calculates likelihood for a sad data set (d=0) assuming correlated errors.
  likelihood                      = ZERO;

  unsigned flagged(0), d(cd), nstop(0);
  double fomreso((double)xtal.sf[d].hires);
  double stopfom(ZERO), meansdluz(ZERO), sdsdluz(ZERO);
  
  if (check)
    xtal.Setselected( (outputmtz ? "PHASE" : "REFINE" ) );

  CMtz::MTZ *MTZOUT(NULL),  *MTZIN(NULL);
  
  if (outputmtz)
  {
    if (allin)
    {
      if (getenv("HKLIN")        != NULL)
        MTZIN                     = CMtz::MtzGet("HKLIN",0);
      else if (xtal.sf[0].mtzin.size())
        MTZIN                     = CMtz::MtzGet(xtal.sf[0].mtzin.c_str(),0);

      if (!MTZIN)
        Bp3Error("Bp3likelihood::sadgradient","HKLIN could not be opened for reading");
    }

    MTZOUT                        = CMtz::MtzMalloc(0,0);

    if (!MTZOUT)
	Bp3Error("Bp3likelihood::sadgradient","HKLOUT could not be opened for writing");
    
    setupmtz(MTZOUT,MTZIN);
  }
  
  unsigned inc                    = (allin) ? colin.size() : 3;

  vector<double> detmodel(xtal.sf[d].nbins, ONE);

  Matrix recovmodel(2), recovinvmodel(2);

  for (unsigned s                 = 0; s < xtal.sf[d].nbins; s++)
  {
    dLdsigmah[d][s]               = dLdsumfpp[d][s]   = ZERO;
    dLdsdluz[d][s]                = ZERO;
    afom[s]                       = cfom[s]           = ZERO;
    anshl[d][s]                   = cnshl[d][s]       = ZERO;
    
    // covariance matrix for the model

    const double model[ ]         = {xtal.sf[d].sigmah[s],
				     xtal.sf[d].sigmah[s]-sumfpp[d][d][s],
				     xtal.sf[d].sigmah[s]-sumfpp[d][d][s],
				     xtal.sf[d].sigmah[s]};
    
    covmodel[s].assign(model);
    inverse2by2(covmodel[s], covinvmodel[s], detmodel[s]);
  }


  for (unsigned r                 = 0; r < xtal.maxselref; r++)
  {
    // Default value for columns to be written out is MNF
    vector<float> fdata(11+inc, CCP4::ccp4_nan().f);
    if (outputhcalc)
      fdata.resize(14+inc, CCP4::ccp4_nan().f);
    
    if (outputmtz)
      storeinitialcolumns(MTZIN, fdata, inc, r);
    
    dLdAp[d][r]                   = dLdBp[d][r] = dLdAm[d][r] =  dLdBm[d][r] = ZERO;

    double eps((double) xtal.epsilon[r]);
    unsigned sa(xtal.bin(d,r));

    // epsilon correct model covariance
    const double tmodel[ ]        = {covmodel[sa](0,0)*eps, covmodel[sa](0,1)*eps,
				     covmodel[sa](1,0)*eps, covmodel[sa](1,1)*eps};

    recovmodel.assign(tmodel);
    
    const double tinvmodel[ ]     = {covinvmodel[sa](0,0)/eps, covinvmodel[sa](0,1)/eps,
				     covinvmodel[sa](1,0)/eps, covinvmodel[sa](1,1)/eps};

    recovinvmodel.assign(tinvmodel);

    double det2(detmodel[sa]*eps*eps);

    double sigmah(xtal.sf[d].sigmah[sa]), sfpp(sumfpp[d][d][sa]);
    double sigman(xtal.sf[d].sigman[sa]), sd(xtal.sf[d].sdluz[sa]);
    unsigned sb(sa);
    double wa(ONE), wb(ZERO);
    
    if (interpolate)
    {
      xtal.binweights(d, r, sa, sb, wa, wb);
      sd                          = (wa*xtal.sf[d].sdluz[sa]    +
				     wb*xtal.sf[d].sdluz[sb]);
    }

    unsigned counts(0);

    if (xtal.sf[d].use(r))
    {
      counts++;
      if (xtal.sf[d].anouse(r) && (outputmtz || !xtal.centric[r]))
	counts++;
    }

    if (counts                    > 1)
    {
      // total covariance matrices

      const double cov[ ]         = {eps*sigman+xtal.sf[d].devp[r]*xtal.sf[d].devp[r],
				     eps*(sigman-sfpp), sd*recovmodel(0,0),
				     sd*recovmodel(0,1),
				     eps*(sigman-sfpp),
				     eps*sigman+xtal.sf[d].devm[r]*xtal.sf[d].devm[r],
				     sd*recovmodel(0,1), sd*recovmodel(0,0),
				     sd*recovmodel(0,0), sd*recovmodel(0,1),
				     recovmodel(0,0), recovmodel(0,1),
				     sd*recovmodel(0,1), sd*recovmodel(0,0),
				     recovmodel(0,1), recovmodel(0,0)};      

      recov.assign(cov);

      double det4(ONE);

      bool filter                 = inverse(recov,recovinv,det4);


      //if (fabs(recov(2,2))        > DSMALL)
      //	sd                        = xtal.sf[d].sdluz[sa] = recov(0,2)/recov(2,2);

      // derivatives of the matrices

      // wrt sdluz

      const double rtmp1[ ]       = {ZERO, ZERO, -recovmodel(0,0), -recovmodel(0,1),
				     ZERO, ZERO, -recovmodel(0,1), -recovmodel(0,0),
				     -recovmodel(0,0), -recovmodel(0,1), ZERO, ZERO,
				     -recovmodel(0,1), -recovmodel(0,0), ZERO, ZERO};
      
      matrixprod(redAdsd,recovinv,rtmp1);

      const double rtmp2[ ]       = {ZERO, ZERO, -eps*sd, -eps*sd,
				     ZERO, ZERO, -eps*sd, -eps*sd,
				     -eps*sd, -eps*sd, -eps, -eps,
				     -eps*sd, -eps*sd, -eps, -eps};

      matrixprod(redAdsigh,recovinv,rtmp2);

      // wrt sumfpp
      const double rtmp3[ ]       = {ZERO, TWO*eps, ZERO, TWO*eps*sd,
				     TWO*eps, ZERO, TWO*eps*sd, ZERO,
				     ZERO, TWO*eps*sd, ZERO, TWO*eps,
				     TWO*eps*sd, ZERO, TWO*eps, ZERO};

      matrixprod(redAdsfpp,recovinv,rtmp3);

      double cosdiff              = tab.Cos(xtal.sf[d].pcalcp[r] - xtal.sf[d].pcalcm[r]);

      double sindiff              = tab.Sin_charged(xtal.sf[d].pcalcp[r] - xtal.sf[d].pcalcm[r]);  

      likelihood                 += (recovinv(0,0)*xtal.sf[d].datap[r]*
				     xtal.sf[d].datap[r]+
				     recovinv(1,1)*xtal.sf[d].datam[r]*
				     xtal.sf[d].datam[r]+
				     (recovinv(2,2) - recovinvmodel(0,0))*
				     xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcp[r] +
				     (recovinv(3,3) - recovinvmodel(1,1))*
				     xtal.sf[d].fcalcm[r]*xtal.sf[d].fcalcm[r]+
				     TWO*xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*
				     (recovinv(2,3) - recovinvmodel(0,1))*cosdiff);

      double cospcalcp            = tab.Cos(xtal.sf[d].pcalcp[r]);
      double sinpcalcp            = tab.Sin_charged(xtal.sf[d].pcalcp[r]);
      double cospcalcm            = tab.Cos(xtal.sf[d].pcalcm[r]);
      double sinpcalcm            = tab.Sin_charged(xtal.sf[d].pcalcm[r]);

      double dldfp                = (TWO*(xtal.sf[d].fcalcp[r]*
					  (recovinv(2,2)-recovinvmodel(0,0)) +
					  xtal.sf[d].fcalcm[r]*
					  (recovinv(2,3)-recovinvmodel(0,1))*
					  cosdiff));
      double dldfm                = (TWO*(xtal.sf[d].fcalcm[r]*
					  (recovinv(3,3)-recovinvmodel(1,1)) +
					  xtal.sf[d].fcalcp[r]*
					  (recovinv(2,3)-recovinvmodel(0,1))*
					  cosdiff));
      double dldpp                = (-TWO*xtal.sf[d].fcalcm[r]*
				     (recovinv(2,3)-recovinvmodel(0,1))*sindiff);
      double dldpm                = (TWO*xtal.sf[d].fcalcp[r]*
				     (recovinv(2,3)-recovinvmodel(0,1))*sindiff);

      
      dLdAp[d][r]                 = dldfp*cospcalcp - dldpp*sinpcalcp;
      dLdBp[d][r]                 = dldfp*sinpcalcp + dldpp*cospcalcp;
      dLdAm[d][r]                 = dldfm*cospcalcm - dldpm*sinpcalcm;
      dLdBm[d][r]                 = dldfm*sinpcalcm + dldpm*cospcalcm;

      double dLdsd                = (redAdsd(0,0)*xtal.sf[d].datap[r]*
				     xtal.sf[d].datap[r]+
				     redAdsd(1,1)*xtal.sf[d].datam[r]*
				     xtal.sf[d].datam[r]+
				     redAdsd(2,2)*xtal.sf[d].fcalcp[r]*
				     xtal.sf[d].fcalcp[r]+
				     redAdsd(3,3)*xtal.sf[d].fcalcm[r]*
				     xtal.sf[d].fcalcm[r]+
				     TWO*xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*
				     redAdsd(2,3)*cosdiff);
      dLdsdluz[d][sa]            += wa*dLdsd;
      
      if (interpolate)
	dLdsdluz[d][sb]          += wb*dLdsd;

      double temp                 = recovinvmodel(0,0) + recovinvmodel(0,1);
      double temp1                = -eps*temp*temp;
      double temp2;
      
      if (updatesigmah)
      {
	dLdsigmah[d][sa]         += (redAdsigh(0,0)*xtal.sf[d].datap[r]*
				     xtal.sf[d].datap[r]+
				     redAdsigh(1,1)*xtal.sf[d].datam[r]*
				     xtal.sf[d].datam[r]+
				     (redAdsigh(2,2) - temp1)*
				     xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcp[r]+
				     (redAdsigh(3,3) - temp1)*
				     xtal.sf[d].fcalcm[r]*xtal.sf[d].fcalcm[r]+
				     TWO*xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*
				     (redAdsigh(2,3) - temp1)*cosdiff);
	
	temp1                     = TWO*eps*(recovinvmodel(0,0)*recovinvmodel(0,0) +
					     recovinvmodel(0,1)*recovinvmodel(0,1));
	temp2                     = FOUR*eps*recovinvmodel(0,0)*recovinvmodel(0,1);

	dLdsumfpp[d][sa]         += (redAdsfpp(0,0)*xtal.sf[d].datap[r]*
				     xtal.sf[d].datap[r]+
				     redAdsfpp(1,1)*xtal.sf[d].datam[r]*
				     xtal.sf[d].datam[r]+
				     (redAdsfpp(2,2) - temp2)*
				     xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcp[r]+
				     (redAdsfpp(3,3) - temp2)*
				     xtal.sf[d].fcalcm[r]*xtal.sf[d].fcalcm[r]+
				     TWO*xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*
				     (redAdsfpp(2,3) - temp1)*cosdiff);
      }
      
      double integral(ZERO), totcos(ZERO), totsin(ZERO); 
      double totlogcos(ZERO), totlogsin(ZERO), totlogcos2(ZERO), totlogsin2(ZERO);
      double dldfcalcp(ZERO), dldfcalcm(ZERO), dldpcalcp(ZERO), dldpcalcm(ZERO);
      double dldsd(ZERO), dldsigh(ZERO), dldsfpp(ZERO);
      double maxexpval(ZERO);

      if (!filter)
      {
	// precalculate arguments in summation/integration      
	double dcos1dfp           = -TWO*xtal.sf[d].datam[r]*cospcalcp*recovinv(1,2);
	double dcos1dfm           = -TWO*xtal.sf[d].datam[r]*cospcalcm*recovinv(1,3);
	double cos1               = (dcos1dfp*xtal.sf[d].fcalcp[r] +
				     dcos1dfm*xtal.sf[d].fcalcm[r]);
	double dcos1dpp           = (TWO*xtal.sf[d].datam[r]*
				     sinpcalcp*recovinv(1,2));
	double dcos1dpm           = (TWO*xtal.sf[d].datam[r]*
				     sinpcalcm*recovinv(1,3)); 
	double dcos1dsd           = (-TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      cospcalcp*redAdsd(1,2) +
				      xtal.sf[d].fcalcm[r]*
				      cospcalcm*redAdsd(1,3)));
	double dcos1dsigh         = (-TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      cospcalcp*redAdsigh(1,2) +
				      xtal.sf[d].fcalcm[r]*
				      cospcalcm*redAdsigh(1,3)));
	double dcos1dsfpp         = (-TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      cospcalcp*redAdsfpp(1,2) +
				      xtal.sf[d].fcalcm[r]*
				      cospcalcm*redAdsfpp(1,3)));

	double dsin1dfp           = (-TWO*xtal.sf[d].datam[r]*
				     sinpcalcp*recovinv(1,2));
	double dsin1dfm           = (-TWO*xtal.sf[d].datam[r]*
				     sinpcalcm*recovinv(1,3));
	double sin1               = (dsin1dfp*xtal.sf[d].fcalcp[r] +
				     dsin1dfm*xtal.sf[d].fcalcm[r]);
	double dsin1dpp           = (-TWO*xtal.sf[d].datam[r]*
				     cospcalcp*recovinv(1,2));
	double dsin1dpm           = (-TWO*xtal.sf[d].datam[r]*
				     cospcalcm*recovinv(1,3));
	double dsin1dsd           = (-TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      sinpcalcp*redAdsd(1,2) +
				      xtal.sf[d].fcalcm[r]*
				      sinpcalcm*redAdsd(1,3)));
	double dsin1dsigh         = (-TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      sinpcalcp*redAdsigh(1,2) +
				      xtal.sf[d].fcalcm[r]*
				      sinpcalcm*redAdsigh(1,3)));
	double dsin1dsfpp         = (-TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      sinpcalcp*redAdsfpp(1,2) +
				      xtal.sf[d].fcalcm[r]*
				      sinpcalcm*redAdsfpp(1,3)));

	temp                      = (recovinv(0,2)*recovinv(0,3)*cosdiff);
      
	double bessarg            = (recovinv(0,1)*recovinv(0,1)*
				     xtal.sf[d].datam[r]*xtal.sf[d].datam[r] +
				     recovinv(0,2)*recovinv(0,2)*
				     xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcp[r] +
				     recovinv(0,3)*recovinv(0,3)*
				     xtal.sf[d].fcalcm[r]*xtal.sf[d].fcalcm[r] +
				     TWO*xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*temp);
      
	double dbessargdfp        = TWO*(xtal.sf[d].fcalcp[r]*
					 recovinv(0,2)*recovinv(0,2) +
					 xtal.sf[d].fcalcm[r]*temp);

	double dbessargdfm        = TWO*(xtal.sf[d].fcalcm[r]*
					 recovinv(0,3)*recovinv(0,3) +
					 xtal.sf[d].fcalcp[r]*temp);

	double dbessargdpp        = -TWO*(xtal.sf[d].fcalcm[r]*
					  recovinv(0,2)*recovinv(0,3)*sindiff);
      
	double dbessargdpm        = TWO*(xtal.sf[d].fcalcp[r]*
					 recovinv(0,2)*recovinv(0,3)*sindiff);
      
	double dbessargdsd        = TWO*(recovinv(0,1)*redAdsd(0,1)*
					 xtal.sf[d].datam[r]*xtal.sf[d].datam[r] +
					 recovinv(0,2)*redAdsd(0,2)*
					 xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcp[r] +
					 recovinv(0,3)*redAdsd(0,3)*
					 xtal.sf[d].fcalcm[r]*xtal.sf[d].fcalcm[r] +
					 xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*
					 ((recovinv(0,2)*redAdsd(0,3) +
					   redAdsd(0,2)*recovinv(0,3))*cosdiff));

	double dbessargdsigh      = TWO*(recovinv(0,1)*redAdsigh(0,1)*
 					 xtal.sf[d].datam[r]*xtal.sf[d].datam[r] +
					 recovinv(0,2)*redAdsigh(0,2)*
					 xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcp[r] +
					 recovinv(0,3)*redAdsigh(0,3)*
					 xtal.sf[d].fcalcm[r]*xtal.sf[d].fcalcm[r] +
					 xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*
					 ((recovinv(0,2)*redAdsigh(0,3) +
					   redAdsigh(0,2)*recovinv(0,3))*cosdiff));
      
	double dbessargdsfpp      = TWO*(recovinv(0,1)*redAdsfpp(0,1)*
					 xtal.sf[d].datam[r]*xtal.sf[d].datam[r] +
					 recovinv(0,2)*redAdsfpp(0,2)*
					 xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcp[r] +
					 recovinv(0,3)*redAdsfpp(0,3)*
					 xtal.sf[d].fcalcm[r]*xtal.sf[d].fcalcm[r] +
					 xtal.sf[d].fcalcp[r]*xtal.sf[d].fcalcm[r]*
					 ((recovinv(0,2)*redAdsfpp(0,3) +
					   redAdsfpp(0,2)*recovinv(0,3))*cosdiff));
      
	temp1                     = (TWO*xtal.sf[d].datam[r]*
				     recovinv(0,1)*recovinv(0,2));
	temp2                     = (TWO*xtal.sf[d].datam[r]*
				     recovinv(0,1)*recovinv(0,3));

	double dcos2dfp           = temp1*cospcalcp;
	double dcos2dfm           = temp2*cospcalcm;
	double cos2               = (dcos2dfp*xtal.sf[d].fcalcp[r] +
				     dcos2dfm*xtal.sf[d].fcalcm[r]);
	double dcos2dpp           = -temp1*sinpcalcp;
	double dcos2dpm           = -temp2*sinpcalcm;

	double dcos2dsd           = (TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      ((recovinv(0,1)*redAdsd(0,2) +
					redAdsd(0,1)*recovinv(0,2))*cospcalcp) +
				      xtal.sf[d].fcalcm[r]*
				      ((recovinv(0,1)*redAdsd(0,3) +
					redAdsd(0,1)*recovinv(0,3))*cospcalcm)));
      
	double dcos2dsigh         = (TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      ((recovinv(0,1)*redAdsigh(0,2) +
					redAdsigh(0,1)*recovinv(0,2))*cospcalcp) +
				      xtal.sf[d].fcalcm[r]*
				      ((recovinv(0,1)*redAdsigh(0,3) +
					redAdsigh(0,1)*recovinv(0,3))*cospcalcm)));

	double dcos2dsfpp         = (TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      ((recovinv(0,1)*redAdsfpp(0,2) +
					redAdsfpp(0,1)*recovinv(0,2))*cospcalcp) +
				      xtal.sf[d].fcalcm[r]*
				      ((recovinv(0,1)*redAdsfpp(0,3) +
					redAdsfpp(0,1)*recovinv(0,3))*cospcalcm)));
      
	double dsin2dfp           = temp1*sinpcalcp;
	double dsin2dfm           = temp2*sinpcalcm;
	double sin2               = (dsin2dfp*xtal.sf[d].fcalcp[r] + 
				     dsin2dfm*xtal.sf[d].fcalcm[r]);
	double dsin2dpp           = temp1*cospcalcp;
	double dsin2dpm           = temp2*cospcalcm;

	double dsin2dsd           = (TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      ((recovinv(0,1)*redAdsd(0,2) +
					redAdsd(0,1)*recovinv(0,2))*sinpcalcp) +
				      xtal.sf[d].fcalcm[r]*
				      ((recovinv(0,1)*redAdsd(0,3) +
					redAdsd(0,1)*recovinv(0,3))*sinpcalcm)));

	double dsin2dsigh         = (TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      ((recovinv(0,1)*redAdsigh(0,2) +
					redAdsigh(0,1)*recovinv(0,2))*sinpcalcp) +
				      xtal.sf[d].fcalcm[r]*
				      ((recovinv(0,1)*redAdsigh(0,3) +
					redAdsigh(0,1)*recovinv(0,3))*sinpcalcm)));

	double dsin2dsfpp         = (TWO*xtal.sf[d].datam[r]*
				     (xtal.sf[d].fcalcp[r]*
				      ((recovinv(0,1)*redAdsfpp(0,2) +
					redAdsfpp(0,1)*recovinv(0,2))*sinpcalcp) +
				      xtal.sf[d].fcalcm[r]*
				      ((recovinv(0,1)*redAdsfpp(0,3) +
					redAdsfpp(0,1)*recovinv(0,3))*sinpcalcm)));

	
	// For stability in the numerical integral, calculate the maximum value
	// of the exponential.  And, since all the values are calculated, store
	// the arguments of exponential and the bessel functions in arrays.

	vector<double> exparg(sadweight.size(), ZERO);
	vector<double> besselarg(sadweight.size(), ZERO);
      
	for (unsigned i           = 0; i < sadweight.size(); i++)
	{ 
	  exparg[i]               = cos1*sadcos[i] + sin1*sadsin[i];
	  double tempbess         = sqrt(std::max(cos2*sadcos[i] +
				 		  sin2*sadsin[i] +
						  bessarg, ZERO));
	  besselarg[i]            = TWO*xtal.sf[d].datap[r]*tempbess;
	  exparg[i]              += besselarg[i];
	  maxexpval               = std::max(maxexpval, exparg[i]);
	}
	
	for (unsigned i           = 0; i < sadweight.size(); i++)
	{
	  double io               = tab.I0e(besselarg[i]);
	  double simarg           = (TWO*xtal.sf[d].datap[r]*xtal.sf[d].datap[r]*
				     tab.Simoverx_preI0(besselarg[i], io));
	  double prob             = tab.ExpM(maxexpval - exparg[i])*io;
	  double wprob            = prob*sadweight[i];
	  integral               += wprob;
	  totcos                 += wprob*sadcos[i];
	  totsin                 += wprob*sadsin[i];

	  dldfcalcp              -= ((dbessargdfp + dcos2dfp*sadcos[i] +
			 	      dsin2dfp*sadsin[i])*simarg +
				     (dcos1dfp*sadcos[i] + dsin1dfp*sadsin[i]))*wprob;
	  dldfcalcm              -= ((dbessargdfm + dcos2dfm*sadcos[i] +
			 	      dsin2dfm*sadsin[i])*simarg +
				     (dcos1dfm*sadcos[i] + dsin1dfm*sadsin[i]))*wprob;
	  dldpcalcp              -= ((dbessargdpp +  dcos2dpp*sadcos[i] +
				      dsin2dpp*sadsin[i])*simarg +
				     (dcos1dpp*sadcos[i] + dsin1dpp*sadsin[i]))*wprob; 
	  dldpcalcm              -= ((dbessargdpm + dcos2dpm*sadcos[i] +
				      dsin2dpm*sadsin[i])*simarg +
				     (dcos1dpm*sadcos[i] + dsin1dpm*sadsin[i]))*wprob;
	  dldsd                  -= ((dbessargdsd + dcos2dsd*sadcos[i] +
			 	      dsin2dsd*sadsin[i])*simarg +
				     (dcos1dsd*sadcos[i] + dsin1dsd*sadsin[i]))*wprob;
	  dldsigh                -= ((dbessargdsigh + dcos2dsigh*sadcos[i] +
			 	      dsin2dsigh*sadsin[i])*simarg +
				     (dcos1dsigh*sadcos[i] + dsin1dsigh*sadsin[i]))*wprob;
	  dldsfpp                -= ((dbessargdsfpp + dcos2dsfpp*sadcos[i] +
				      dsin2dsfpp*sadsin[i])*simarg +
				     (dcos1dsfpp*sadcos[i] + dsin1dsfpp*sadsin[i]))*wprob;
	  if (prob                > SMALLESTD)
	  {
	    prob                  = log(prob)*sadweight[i];
	    totlogcos            += prob*sadcos[i];
	    totlogsin            += prob*sadsin[i];
	    totlogcos2           += prob*(sadcos[i]*sadcos[i] - sadsin[i]*sadsin[i]);
	    totlogsin2           += prob*TWO*sadcos[i]*sadsin[i];
	  }
	}
      }
      
      if ( (integral              > DSMALL) && !filter) 
      {
	likelihood               -= (log(det2*integral/det4) + maxexpval);

	dLdAp[d][r]              += (dldfcalcp*cospcalcp - 
				     dldpcalcp*sinpcalcp)/integral;
	dLdBp[d][r]              += (dldfcalcp*sinpcalcp + 
				     dldpcalcp*cospcalcp)/integral;
	dLdAm[d][r]              += (dldfcalcm*cospcalcm - 
				     dldpcalcm*sinpcalcm)/integral;
	dLdBm[d][r]              += (dldfcalcm*sinpcalcm + 
				     dldpcalcm*cospcalcm)/integral;

       double dLdsd               = (dldsd/integral +
				     TWO*((recovmodel(0,0)*(recovinv(0,2) +
							    recovinv(1,3)) +
					   recovmodel(0,1)*(recovinv(0,3) +
							    recovinv(1,2)))));
	
       dLdsdluz[d][sa]          += wa*dLdsd;

       if (interpolate)
	 dLdsdluz[d][sb]        += wb*dLdsd;

       if (updatesigmah)
       {
	 dLdsigmah[d][sa]       += (dldsigh/integral +
				    eps*((recovinv(0,2) + recovinv(0,3) +
	 				   recovinv(1,2) + recovinv(1,3))*TWO*sd +
					  recovinv(2,2) + TWO*recovinv(2,3) +
					  recovinv(3,3) -
					  (recovinvmodel(0,0)+TWO*recovinvmodel(0,1) +
					   recovinvmodel(1,1))));
	
	 dLdsumfpp[d][sa]       += (dldsfpp/integral -
				     FOUR*eps*(recovinv(0,1) + recovinv(2,3) +
					       (recovinv(0,3) + recovinv(1,2))*sd -
					       recovinvmodel(0,1)));	
       }
       
       double fom                 = sqrt(totcos*totcos + totsin*totsin)/integral;

	if (xtal.Getres(d,r)     >= fomreso)
	{
	  stopfom                += fom;
	  meansdluz              += sd;
	  sdsdluz                += sd*sd;
	  nstop++;
	}

	double phib               = atan2(totsin, totcos);
	if (xtal.centric[r])
	  phib                    =  ((cos(phib - xtal.centricphase[r]) >= ZERO) ?
				      xtal.centricphase[r] : xtal.centricphase[r] + PI);
	  
	unsigned s                = xtal.bin(d,r);

	if (xtal.centric[r])
	{
	  cfom[s]                += fom;
	  cnshl[d][s]++;
	}
	else
	{
	  afom[s]                += fom;
	  anshl[d][s]++;
	}
	  
	if (outputmtz)
	{
	  fdata[inc + 2]          = (float) (fom*xtal.sf[d].datap[r]);
	  fdata[inc + 3]          = (float) (phib/DEGREEtoRAD);
	  fdata[inc + 4]          = (float) fom;
	  fdata[inc + 5]          = (float) totlogcos;
	  fdata[inc + 6]          = (float) totlogsin;
	  fdata[inc + 7]          = (float) totlogcos2;
	  fdata[inc + 8]          = (float) totlogsin2;
	  double fpp(mdl.form[mdl.atom[0].nform].fpp[0]);
	  double adiff            =  fpp*(-dLdBp[d][r] + dLdBm[d][r]);
	  double bdiff            =  fpp*( dLdAp[d][r] - dLdAm[d][r]);
	  fdata[inc + 9]          = (float) sqrt(adiff*adiff + bdiff*bdiff);
	  fdata[inc + 10]         = (float) (atan2(bdiff,adiff)/DEGREEtoRAD);
	}
      }
      else
      {
	// double tmp(std::max(fabs(det4-0.001), ONE));
	// likelihood                += 1000.0*(log(tmp) + 1.75);
	likelihood                += 2750.0;
/*
	if (tmp                    > ONE)
	{
	  dLdsdluz[d][sa]         += 1000.0*det4*(TWO*((recovinv(0,2)+recovinv(1,3))*recov(2,2) +
                                                       (recovinv(0,3)+recovinv(1,2))*recov(2,3)))/(0.001-det4);
          dLdsigmah[d][sa]        += 1000.0*det4*eps*(TWO*sd*(recovinv(0,2)+recovinv(0,3)+recovinv(1,2)+
                                                             recovinv(1,3))+recovinv(2,2)+recovinv(3,3)+
                                                      TWO*recovinv(2,3))/(0.001-det4);
          dLdsumfpp[d][sa]        += 1000.0*det4*(-FOUR*eps*(recovinv(0,1)+sd*(recovinv(0,3)+recovinv(1,2))+
                                                             recovinv(2,3)))/(0.001-det4);
        }
*/	
	flagged++;
      }
    }
    else if (outputmtz && (xtal.sf[d].usem(r) || xtal.sf[d].usep(r)) )
    {
      double fobs, arg(ZERO);
      double var                  = std::max(eps*(sigman - sd*sd*sigmah), EPSILON);
	
      if (xtal.sf[d].usep(r))
      {
	fobs                      = xtal.sf[d].datap[r];
	var                      += xtal.sf[d].devp[r]*xtal.sf[d].devp[r];
	arg                       = TWO*fobs*sd*xtal.sf[d].fcalcp[r]/var;
      }
      else
      {
	fobs                      = xtal.sf[d].datam[r];
	var                      += xtal.sf[d].devm[r]*xtal.sf[d].devm[r];
	arg                       = TWO*fobs*sd*xtal.sf[d].fcalcm[r]/var;
      }
	
      double phib                 = xtal.sf[d].pcalcp[r];

      double fom                  = tab.Sim(arg);

      unsigned s                  = xtal.bin(d,r);
      afom[s]                    += fom;
      anshl[d][s]++;
 
      fdata[inc + 2]              = (float) fom*fobs;
      fdata[inc + 3]              = (float) (phib/DEGREEtoRAD);
      fdata[inc + 4]              = (float) fom;
      double temp                 = atanh(fom);
      fdata[inc + 5]              = (float) temp*tab.Cos(phib);
      fdata[inc + 6]              = (float) temp*tab.Sin_charged(phib);
      fdata[inc + 7]              = (float) ZERO;
      fdata[inc + 8]              = (float) ZERO;
      fdata[inc + 9]              = (float) ZERO;
      fdata[inc + 10]             = (float) ZERO;
    }
    
    if (outputmtz)
    {      
      if (outputhcalc)
      {
	unsigned sa(xtal.bin(d,r)), sb(sa);
	double wa(ONE), wb(ZERO);
	double sig               = wa*xtal.sf[d].sigmah[sa] + wb*xtal.sf[d].sigmah[sb];
	xtal.binweights(d, r, sa, sb, wa, wb);

	fdata[inc + 11]              = (float) xtal.sf[d].fcalcp[r];
	fdata[inc + 12]              = (float) xtal.sf[d].pcalcp[r]/DEGREEtoRAD;
	fdata[inc + 13]              = (float) xtal.sf[d].fcalcp[r]/sqrt(eps*sig);
      }
      CMtz::ccp4_lwrefl(MTZOUT, &fdata[0], &colout[0], fdata.size(), r+1);
    }
  }

  if (outputmtz)
  {
    if (nstop)
    {
      double dnstop((double)nstop);
      char result[125];
      sprintf(result,"The overall FOM is %.3f and the average anomalous Luzzati error is %.3f to %.2f resolution",
 	      stopfom/dnstop,meansdluz/dnstop,fomreso);
      Bp3Result(result);
      printf("The standard deviation of the anomalous Luzzati error is %.3f to %.2f resolution\n\n",
	     sqrt(fabs(dnstop*sdsdluz - meansdluz*meansdluz)/(dnstop*(dnstop-1))), fomreso);
      FILE *pluz(NULL);
      pluz = fopen("pluzz1","w");
      fprintf(pluz," PLUZzati 2" );
      for (unsigned s = 0; s < xtal.sf[d].nbins; s++)
	fprintf(pluz, "%9.6f ",xtal.sf[d].sdluz[s]);
      fprintf(pluz, "NOREF \n");
//      fprintf(pluz, " \n");
      fclose(pluz);

    }
    
    CMtz::MtzPut(MTZOUT, " ");
    if (MTZOUT)
      CMtz::MtzFree(MTZOUT);
    if (MTZIN)
      CMtz::MtzFree(MTZIN);

  }  
  
  if (xtal.verbose                > 1)
    printf("LIKELIHOOD = %f\n",likelihood);
  
  if (flagged && check)
    printf("Number of reflections with low likelihood values: %u\n", flagged);

  Tsadgradient.stop("Main");
  return likelihood;
}
