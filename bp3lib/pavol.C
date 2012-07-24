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
#include "likelihood.h"
#include <complex>
#include <fstream>
#include <iostream>
#include "llhood.h"
#include "covmat.h"
#include "gaussian.h"
#include <time.h>

typedef double Type;

double Likelihood::multsirasgradient(const bool check, const bool outputmtz)
{
  likelihood                      = ZERO;
  bool phasecomb((protocol       == "PHASECOMB"));
  
// CONSTRUCTION OF THE LIKELIHOOD CLASS - this should be done only once at the beginning of the run!
  int Num(5);
  int N_meas(3);
  int N_part(xtal.pluz.size()+2);
  /*
  if (updatesigmah)
    N_part                       += 2;
  */ 
  
 // the extra 2 are not real D's but only used to get derivs wrt sumfpp and sigmah
  int N_mod(Num - N_meas);
  string targe("srasph");
  multivar_llhood::likelihood<Type> l( N_meas, Num, 1, 1, targe.c_str(), N_part, 100, 1000, 1 );
  l.checkOK                       = 1;
  l.SetImproveInteg(4); 
  l.cov->no_imag                  = 1;
  l.SetNoIntegLast(0);

  Type *F                         = new Type[Num];
  Type *ph                        = new Type[Num];

  unsigned flagged(0), nstop(0);
  double fomreso((double)xtal.sf[0].hires);
  double stopfom(ZERO);
  
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
	Bp3Error("Likelihood::multsirasgradient","HKLIN could not be opened for reading");
    }

    MTZOUT                        = CMtz::MtzMalloc(0,0);

    if (!MTZOUT)
      Bp3Error("Likelihood::multsirasgradient","HKLOUT could not be opened for writing");
    
    setupmtz(MTZOUT,MTZIN);
    
  }
  
  unsigned inc                    = (allin) ? colin.size() : 3;

  for (unsigned p                 = 0; p < xtal.pluz.size(); p++)
    for (unsigned s               = 0; s < xtal.sf[1].nbins; s++)
      dLdpluz[p][s]               = ZERO;
  
  for (unsigned d                 = 0; d < xtal.sf.size()  ; d++)
  {
    //    initstat();
    for (unsigned s               = 0; s < xtal.sf[d].nbins; s++)
    {
      dLdsigmah[d][s]             = dLdsumfpp[d][s]   = ZERO;
      afom[s]                     = cfom[s]           = ZERO;
      anshl[d][s]                 = cnshl[d][s]       = ZERO;
    }  
  }
  
  
  // start the evaluations   	
  int neg_eig_num = 0;
  for (unsigned r                 = 0; r < xtal.maxselref; r++)
  {
    // Default value for columns to be written out is MNF
    vector<float> fdata(9+inc, CCP4::ccp4_nan().f);
    
    if (outputmtz)
  	  storeinitialcolumns(MTZIN, fdata, inc, r);
    
    double eps((double) xtal.epsilon[r]);
    unsigned sa1(xtal.bin(1,r));
  	
    unsigned counts(0);
      
    for (unsigned d               = 0; d < xtal.sf.size(); d++)
    {
      dLdAp[d][r]                 = dLdBp[d][r] = ZERO;
      if (xtal.sf[d].anomalous)
	dLdAm[d][r]               = dLdBm[d][r] = ZERO;
  	
      if (xtal.sf[d].use(r))
      {
	counts++;
	if (xtal.sf[d].anomalous)
	  if (xtal.sf[d].anouse(r) || xtal.centric[r])
	    counts++;
      }
    }
      
    if  ( counts                  > 2)
    {
      F[0]                        = xtal.sf[0].datap[r];
      F[1]                        = xtal.sf[1].datap[r];
      F[2]                        = xtal.sf[1].datam[r];
      Type Ap(xtal.sf[1].fcalcp[r]*cos(xtal.sf[1].pcalcp[r]));
      Type Bp(xtal.sf[1].fcalcp[r]*sin(xtal.sf[1].pcalcp[r]));
      Type Am(xtal.sf[1].fcalcm[r]*cos(xtal.sf[1].pcalcm[r]));
      Type Bm(xtal.sf[1].fcalcm[r]*sin(xtal.sf[1].pcalcm[r]));
      Type A3((Ap+Am)/2);
      Type B3((Bp+Bm)/2);
      Type A4((Ap-Am)/2);
      Type B4((Bp-Bm)/2);
      F[3]                        = sqrt(A3*A3+B3*B3);
      F[4]                        = sqrt(A4*A4+B4*B4);
      ph[3]                       = atan2(B3,A3);
      ph[4]                       = atan2(B4,A4);
  	
      // INPUT FOR SIRAS
  	
      l.SetCentRice(0,0);
  	
//   measured sigmas
      l.cov->sigma_N[0]           = xtal.sf[1].sigmanref[sa1]*eps;
      l.cov->sigma_N[1]           = xtal.sf[1].sigman[sa1]*eps;
      l.cov->sigma_N[2]           = (xtal.sf[1].sigman[sa1] - sumfpp[1][1][sa1])*eps;
  	
//   calculated sigmas
      l.cov->part[0].sigma_P[0]   = 0.;
      l.cov->part[0].sigma_P[1]   = xtal.sf[1].sigmah[sa1]*eps;
      l.cov->part[0].sigma_P[2]   = sumfpp[1][1][sa1]/2.*eps;
        
//   meas errors
      l.cov->sig_meas[0]          = xtal.sf[0].devp[r]; 
      l.cov->sig_meas[1]          = xtal.sf[1].devp[r]; 
      l.cov->sig_meas[2]          = xtal.sf[1].devm[r];

      int phibcalc((int)outputmtz-1);
      l.SetMaxPoints(10000000);
      if (phibcalc               >= 0) l.SetMaxPoints(10000000);
      if (phibcalc               >= 0) l.SetImproveInteg(3);
      if (phibcalc               >= 0) l.SetSaferNumPoints(3);
      l.SetNum4PHIB(phibcalc);
      l.SetCalcHL(phibcalc+1);

// D's
      for (unsigned p             = 0; p < xtal.pluz.size(); p++)
	l.cov->part[p].D[0][1]    = l.cov->part[p].D[1][0] = xtal.pluz[p][sa1];
	

      // CALL THE FUNCTIONS TO CALCULATE MATRIX ANC ITS INVERSE AND EIGENVALUES
      l.cov->Make_matrix();
      l.InverseAndEigen();

      int out(0);
      for (int i                  = 0; i < l.Num; i++)
	if ( l.eigenvalues1[l.Rice][i]   < l.GetMinEig() )
	  out                     = 1;
      if (!out)
	likelihood               += l.EvaluateSR( F, ph );
      else
      {
    neg_eig_num ++;
	likelihood               += l.EvaluateSR( F, ph ) + 2500.0;
// only luzzati commented out to resolve non-positive definiteness problems in !onlyluzzati cases
	if (check)
	{
 	  for (unsigned s         = 0; s < xtal.sf[1].nbins; s++)
 	  {
	    if (xtal.refinep[0])
	      xtal.pluz[0][s]  *= 0.95;
	    if (xtal.refinep[1])
	      xtal.pluz[1][s]  *= 0.995;
	    if (xtal.refinep[2])
	      xtal.pluz[2][s]  *= 1.001;
 	    if (xtal.refinep[3])
	      xtal.pluz[3][s]    *= 0.95;
 	  }
	}	
      }
      if (!out)
      {
	// transforming derivs wrt F,ph to derivs wrt A,B 
	int Nm3(3);
	Type F_Nm3(F[Nm3]);
	Type dF_dA3(0.), dF_dB3(0.), dph_dA3(0.), dph_dB3(0.);
	if (F_Nm3                 > 1e-7)
	{
	  dF_dA3                  = A3/F_Nm3;
	  dF_dB3                  = B3/F_Nm3;
	  dph_dA3                 = - B3/F_Nm3/F_Nm3;
	  dph_dB3                 = A3/F_Nm3/F_Nm3;
	}
	Type DFDA3(l.GetDer_F(Nm3)*dF_dA3 + l.GetDer_ph(Nm3)*dph_dA3);
	Type DFDB3(l.GetDer_F(Nm3)*dF_dB3 + l.GetDer_ph(Nm3)*dph_dB3);
	
	int Nm4(Nm3+1);
	Type F_Nm4(F[Nm4]);
	Type dF_dA4(0.), dF_dB4(0.), dph_dA4(0.), dph_dB4(0.);
	if (F_Nm4                 > 1e-7)
	{
	  dF_dA4                  = A4/F_Nm4;
	  dF_dB4                  = B4/F_Nm4;
	  dph_dA4                 = - B4/F_Nm4/F_Nm4;
	  dph_dB4                 = A4/F_Nm4/F_Nm4;
	}
	Type DFDA4(l.GetDer_F(Nm4)*dF_dA4 + l.GetDer_ph(Nm4)*dph_dA4);
	Type DFDB4(l.GetDer_F(Nm4)*dF_dB4 + l.GetDer_ph(Nm4)*dph_dB4);

	//transforming derivs wrt A,B anom, non-anom to derivs wrt A,B +,-
	dLdAp[1][r]               = HALF*(DFDA3+DFDA4);
      
	dLdBp[1][r]               = HALF*(DFDB3+DFDB4);
      
	dLdAm[1][r]               = HALF*(DFDA3-DFDA4);
      
	dLdBm[1][r]               = HALF*(DFDB3-DFDB4);
	for (unsigned p           = 0; p < xtal.pluz.size(); p++)
	  dLdpluz[p][sa1]        += l.GetDer_D(p,0,1);

	if (updatesigmah)
	{ 
	   dLdsigmah[1][sa1]     += l.GetDer_D( xtal.pluz.size(), 0,1)*eps;
      
	   dLdsumfpp[1][sa1]     += l.GetDer_D( xtal.pluz.size()+1, 0,1)*eps;
	}
	
        double fom                = l.GetFOM();

	if (xtal.Getres(0,r)     >= fomreso)
	{
	  stopfom                += fom;
	  nstop++;
	}
      
	double phib               = l.GetPHIB();
	if (xtal.centric[r])
	  phib                    =  ((cos(phib - xtal.centricphase[r]) >= ZERO) ?
				      xtal.centricphase[r] : xtal.centricphase[r] + PI);
	
	if (xtal.centric[r])
	  cfom[xtal.bin(0,r)]    += fom;
	else
	  afom[xtal.bin(0,r)]    += fom;

	/*
	for (unsigned d           = 0; d < xtal.sf.size()          ; d++) 
	  sanonshl[d][xtal.bin(d,r)]++;
	*/

	for (unsigned d           = 0; d < xtal.sf.size()          ; d++) 
	  if (xtal.centric[r])
	    cnshl[d][xtal.bin(d,r)]++;
	  else
	    anshl[d][xtal.bin(d,r)]++;
	  
	if (outputmtz)
	{
	  fdata[inc + 2]          =  (float) (fom*xtal.sf[0].datap[r]);
	  fdata[inc + 3]          = (float) (phib/DEGREEtoRAD);
	  fdata[inc + 4]          = (float) fom;
	  fdata[inc + 5]          = (float) l.GetHLA();  // HLA
	  fdata[inc + 6]          = (float) l.GetHLB();  // HLB
	  fdata[inc + 7]          = (float) l.GetHLC();  // HLC
	  fdata[inc + 8]          = (float) l.GetHLD();  // HLD
     	}
      }
      else { /*out!*/ }
    }
    
    if (outputmtz)
      CMtz::ccp4_lwrefl(MTZOUT, &fdata[0], &colout[0], fdata.size(), r+1);
  }

  if (outputmtz)
  {
    if (nstop)
    {
      double dnstop((double)nstop);
      string result;
      sprintf(&result[0],"The overall FOM is %.3f\n",
 	      stopfom/dnstop);
      Bp3Result(result.c_str());
    }
    
    CMtz::MtzPut(MTZOUT, " ");
    if (MTZOUT)
      CMtz::MtzFree(MTZOUT);
    if (MTZIN)
      CMtz::MtzFree(MTZIN);
  }  
  
  if (xtal.verbose                > 0)
    printf("LIKELIHOOD = %f\n",likelihood);

  if (xtal.verbose                > 0)
    cout << "NUMBER OF NEG EIGENVALUES = " << neg_eig_num << endl;
  
  if (flagged && check)
    printf("Number of reflections with low likelihood values: %u\n", flagged);

  delete[] F;
  delete[] ph;
  
  return likelihood;
}


// SAD PART

double Likelihood::pavolsadgradient(const bool check, const bool outputmtz)
{
  likelihood                      = ZERO;
  bool phasecomb((protocol       == "PHASECOMB"));
  
// CONSTRUCTION OF THE LIKELIHOOD CLASS - this should be done only once at the beginning of the run!
  int Num(4);
  if (phasecomb)
    Num++;

  int N_meas(2);

  int N_part(xtal.pluz.size()+2);

  /*
  if (updatesigmah)
    N_part                       += 2; // the extra 2 are not real D's but only used to get derivs wrt sumfpp and sigmah
  */

  int N_mod(Num - N_meas);
  string targe("sadh");
  
  multivar_llhood::likelihood<Type> l( N_meas, Num, 1, 1, targe.c_str(), N_part, 25, 1, 1 );
  if (target                     == "PSAD")
    l.cov->SetZeroRows(4);
  l.checkOK                       = 1;
  l.SetImproveInteg(0); 
  l.cov->no_imag                  = 1;
  l.SetNoIntegLast(0);

  Type *F                         = new Type[Num];
  Type *ph                        = new Type[Num];

  unsigned flagged(0), nstop(0);
  double fomreso((double)xtal.sf[0].hires);
  double stopfom(ZERO);
  
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
	Bp3Error("Likelihood::pavolgradient","HKLIN could not be opened for reading");
    }

    MTZOUT                        = CMtz::MtzMalloc(0,0);

    if (!MTZOUT)
      Bp3Error("Likelihood::pavolgradient","HKLOUT could not be opened for writing");
    
    setupmtz(MTZOUT,MTZIN);
    
  }
  
  unsigned inc                    = (allin) ? colin.size() : 3;

  for (unsigned p                 = 0; p < xtal.pluz.size(); p++)
    for (unsigned s               = 0; s < xtal.sf[0].nbins; s++)
      dLdpluz[p][s]               = ZERO;
  
  for (unsigned d                 = 0; d < xtal.sf.size()          ; d++)
  {
    // initstat();
    for (unsigned s               = 0; s < xtal.sf[d].nbins; s++)
    {
      dLdsigmah[d][s]             = dLdsumfpp[d][s]   = ZERO;
      afom[s]                     = cfom[s]           = ZERO;
      anshl[d][s]                 = cnshl[d][s]       = 0;
    }  
  }
  
  
  // start the evaluations   	
  for (unsigned r                 = 0; r < xtal.maxselref; r++)
  {
    // Default value for columns to be written out is MNF
    unsigned cols(11);
    if (phasecomb)
      cols                       += 2;
    
    vector<float> fdata(cols+inc, CCP4::ccp4_nan().f);
    
    if (outputmtz)
      storeinitialcolumns(MTZIN, fdata, inc, r);
    
    double eps((double) xtal.epsilon[r]);
    unsigned sa1(xtal.bin(0,r));
  	
    unsigned counts(0);
      
    for (unsigned d               = 0; d < xtal.sf.size(); d++)
    {
      dLdAp[d][r]                 = dLdBp[d][r] = ZERO;
      if (xtal.sf[d].anomalous)
	dLdAm[d][r]               = dLdBm[d][r] = ZERO;
  	
      if (xtal.sf[d].use(r))
      {
	counts++;
	if (xtal.sf[d].anomalous)
	  if (xtal.sf[d].anouse(r) && (outputmtz || !xtal.centric[r]))
	    counts++;
      }
    }
      
    if  ( counts                  > 1)
    {
      F[0]                        = xtal.sf[0].datap[r];
      F[1]                        = xtal.sf[0].datam[r];
      Type Ap                     = xtal.sf[0].fcalcp[r]*cos(xtal.sf[0].pcalcp[r]);
      Type Bp                     = xtal.sf[0].fcalcp[r]*sin(xtal.sf[0].pcalcp[r]);
      Type Am                     = xtal.sf[0].fcalcm[r]*cos(xtal.sf[0].pcalcm[r]);
      Type Bm                     = xtal.sf[0].fcalcm[r]*sin(xtal.sf[0].pcalcm[r]);
      Type A4                     = (Ap+Am)/2;
      Type B4                     = (Bp+Bm)/2;
      Type A3                     = (Ap-Am)/2;
      Type B3                     = (Bp-Bm)/2;
      F[2]                        = sqrt(A3*A3+B3*B3);
      ph[2]                       = atan2(B3,A3);
      l.cov->part[0].sigma_P[1]   = sumfpp[0][0][sa1]/2.*eps;
      if (phasecomb)
      {
	F[3]                      = xtal.sf[0].fmodel[r][0];
	F[4]                      = sqrt(A4*A4+B4*B4);
	ph[3]                     = xtal.sf[0].pmodel[r][0]*DEGREEtoRAD;
	ph[4]                     = atan2(B4,A4);
	l.cov->part[0].sigma_P[0] = xtal.sf[0].sigmap[sa1]*eps; // F[3]
	l.cov->part[0].sigma_P[2] = xtal.sf[0].sigmah[sa1]*eps; // F[4]
	// printf("sumfpp = %f sigmah = %f\n", sumfpp[0][0][sa1]/2.*eps, xtal.sf[0].sigmah[sa1]*eps);
      }
      else
      {
	F[3]                      = sqrt(A4*A4+B4*B4);
	ph[3]                     = atan2(B4,A4);
	l.cov->part[0].sigma_P[0] = xtal.sf[0].sigmah[sa1]*eps;
      }
  	
      // INPUT FOR SAD
  	
      l.SetCentRice(0,0);
  	
//   measured sigmas
      l.cov->sigma_N[0]           = xtal.sf[0].sigman[sa1]*eps;
      l.cov->sigma_N[1]           = (xtal.sf[0].sigman[sa1] - sumfpp[0][0][sa1])*eps;  	
        
//   meas errors
      l.cov->sig_meas[0]          = xtal.sf[0].devp[r]; 
      l.cov->sig_meas[1]          = xtal.sf[0].devm[r];

      int phibcalc                = (int)outputmtz-1;
//      l.SetNum4PHIB(phibcalc);
      l.SetCalcHL(phibcalc+1);

// D's
      for (unsigned p             = 0; p < xtal.pluz.size(); p++)
	    l.cov->part[p].D[0][1]  = l.cov->part[p].D[1][0] = xtal.pluz[p][sa1];

      // CALL THE FUNCTIONS TO CALCULATE MATRIX ANC ITS INVERSE AND EIGENVALUES
      l.cov->Make_matrix();
      l.InverseAndEigen();
      
      int out(0);
      for (int i                  = 0; i < l.Num; i++)
	if ( l.eigenvalues1[l.Rice][i] < l.GetMinEig() )
	  out                     = 1;
      Type FVALUE                 = l.EvaluateSR( F, ph );

      /*
      if (xtal.verbose)
        l.cov->Print();
      */
      
//cout << r << " "<< xtal.miller[r][0] <<" "<<xtal.miller[r][1] << " " << xtal.miller[r][2] << endl;
//cout << F[0] << " " << F[1] << " " << F[2] << " " << F[3] << " " << ph[2]<< " "<< ph[3] << endl;
//cout << l.cov->sigma_N[0] << " "<< l.cov->sigma_N[1] <<" "<< l.cov->part[0].sigma_P[0] << " "<< l.cov->part[0].sigma_P[1] << endl;
//cout <<l.cov->part[0].D[0][1] << " "<< l.cov->part[1].D[0][1] << " "<< l.cov->part[2].D[0][1] << " "<< l.cov->part[3].D[0][1] << " "<< endl;
//cout << FVALUE << " "<<l.GetPHIB() << endl;
      if (!out)
	likelihood               += FVALUE;
      else
      {
	likelihood               += FVALUE + 2500.0;
	flagged++;
	if (check)
	  for (unsigned s         = 0; s < xtal.sf[0].nbins; s++)
	  {
	    xtal.pluz[0][s]      *= 0.95;  // correlation between obs and dm 
	    xtal.pluz[1][s]      *= 0.995; // anomalous luzzati
	    // xtal.pluz[2][s]      *= 1.001; // scale (scale between observation and scale)  ***NSP
	    //	    xtal.pluz[3][s]      *= 0.95;
	  }
      }
      
//cout << FVALUE << endl;

      if (!out)
      {
// transforming derivs wrt F,ph to derivs wrt A,B 
	int Nm3                   = 2;
	Type F_Nm3                = F[Nm3];
	Type dF_dA3(0.), dF_dB3(0.), dph_dA3(0.), dph_dB3(0.);
	if (F_Nm3                 > 1e-7)
	{
	  dF_dA3                  = A3/F_Nm3;
	  dF_dB3                  = B3/F_Nm3;
	  dph_dA3                 = - B3/F_Nm3/F_Nm3;
	  dph_dB3                 = A3/F_Nm3/F_Nm3;
	}
	Type DFDA3                = l.GetDer_F(Nm3)*dF_dA3 + l.GetDer_ph(Nm3)*dph_dA3;
	Type DFDB3                = l.GetDer_F(Nm3)*dF_dB3 + l.GetDer_ph(Nm3)*dph_dB3;
	
	int Nm4                   = Nm3+1;
	Type F_Nm4                = F[Nm4];
	Type dF_dA4(0.), dF_dB4(0.), dph_dA4(0.), dph_dB4(0.);
	if (F_Nm4                 > 1e-7) {
	  dF_dA4                  = A4/F_Nm4;
	  dF_dB4                  = B4/F_Nm4;
	  dph_dA4                 = - B4/F_Nm4/F_Nm4;
	  dph_dB4                 = A4/F_Nm4/F_Nm4;
	}
	Type DFDA4                = l.GetDer_F(Nm4)*dF_dA4 + l.GetDer_ph(Nm4)*dph_dA4;
	Type DFDB4                = l.GetDer_F(Nm4)*dF_dB4 + l.GetDer_ph(Nm4)*dph_dB4;

	//transforming derivs wrt A,B anom, non-anom to derivs wrt A,B +,-
	dLdAp[0][r]               = HALF*(DFDA3+DFDA4);
      
	dLdBp[0][r]               = HALF*(DFDB3+DFDB4);
      
	dLdAm[0][r]               = HALF*(-DFDA3+DFDA4);
      
	dLdBm[0][r]               = HALF*(-DFDB3+DFDB4);
	
	for (unsigned p           = 0; p < xtal.pluz.size(); p++)
	  dLdpluz[p][sa1]        += l.GetDer_D(p,0,1);

	if (updatesigmah)
	{
	  // ***NSP
	  dLdsigmah[0][sa1]      += l.GetDer_D( xtal.pluz.size(), 0,1)*eps;
      
	  dLdsumfpp[0][sa1]      += l.GetDer_D( xtal.pluz.size()+1, 0,1)*eps;
	}
	
        double fom                = l.GetFOM();

	if (xtal.Getres(0,r)     >= fomreso)
	{
	  stopfom                += fom;
	  nstop++;
	}
      
	double phib               = l.GetPHIB();
	if (xtal.centric[r])
	  phib                    =  ((cos(phib - xtal.centricphase[r]) >= ZERO) ?
				      xtal.centricphase[r] : xtal.centricphase[r] + PI);
	
	if (xtal.centric[r])
	{
	  cfom[sa1]              += fom;
	  cnshl[0][sa1]++;
	}
	else
	{
	  afom[sa1]              += fom;
	  anshl[0][sa1]++;
	}	
	  
	if (outputmtz)
	{
	  double fobs             = (xtal.sf[0].datap[r]+xtal.sf[0].datam[r])*HALF;
	  fdata[inc + 2]          = (float) (fom*fobs);
	  fdata[inc + 3]          = (float) (phib/DEGREEtoRAD);
	  fdata[inc + 4]          = (float) fom;
	  fdata[inc + 5]          = (float) l.GetHLA();  // HLA
	  fdata[inc + 6]          = (float) l.GetHLB();  // HLB
	  fdata[inc + 7]          = (float) l.GetHLC();  // HLC
	  fdata[inc + 8]          = (float) l.GetHLD();  // HLD

	  double fpp(mdl.form[mdl.atom[0].nform].fpp[0]);
	  double adiff            =  fpp*(-dLdBp[0][r] + dLdBm[0][r]);
	  double bdiff            =  fpp*( dLdAp[0][r] - dLdAm[0][r]);
	  fdata[inc + 9]          = (float) sqrt(adiff*adiff + bdiff*bdiff);
	  fdata[inc + 10]         = (float) (atan2(bdiff,adiff)/DEGREEtoRAD);
	  if (phasecomb)
 	  {
	    double acomb          = TWO*fom*fobs*cos(phib) - xtal.pluz[0][sa1]*F[3]*cos(ph[3]);
	    double bcomb          = TWO*fom*fobs*sin(phib) - xtal.pluz[0][sa1]*F[3]*sin(ph[3]);
	    // fcomb and phicomb ***NSP
	    fdata[inc + 11]       = (float) (sqrt(acomb*acomb + bcomb*bcomb));
	    fdata[inc + 12]       = (float) (atan2(bcomb,acomb)/DEGREEtoRAD);
	  }
	}
      }
      else { /*out!*/ }
    }
    else if (outputmtz)
    {
      double dluz(xtal.pluz[0][sa1]), fcalc(xtal.sf[0].fmodel[r][0]), pcalc(xtal.sf[0].pmodel[r][0]*DEGREEtoRAD);
      if (xtal.sf[0].use(r))
      {
	double fobs, arg(ZERO), dluz(xtal.pluz[0][sa1]);
	double var                = std::max(eps*(xtal.sf[0].sigman[sa1] - 
						  dluz*dluz*xtal.sf[0].sigmap[sa1]), EPSILON);
	if (xtal.sf[0].usep(r))
	{
	  fobs                    = xtal.sf[0].datap[r];
	  var                    += xtal.sf[0].devp[r]*xtal.sf[0].devp[r];
	  arg                     = TWO*fobs*dluz*fcalc/var;
	}
	else
	{
	  fobs                    = xtal.sf[0].datam[r];
	  var                    += xtal.sf[0].devm[r]*xtal.sf[0].devm[r];
	  arg                     = TWO*fobs*dluz*fcalc/var;
	}

	double fom                = tab.Sim(arg);

	if (xtal.centric[r])
	{
	  cfom[sa1]                += fom;
	  cnshl[0][sa1]++;
	}
	else
	{
	  afom[sa1]                += fom;
	  anshl[0][sa1]++;
	}	

	fdata[inc + 2]          = (float) (fom*fobs);
	fdata[inc + 3]          = (float) (xtal.sf[0].pmodel[r][0]);
	fdata[inc + 4]          = (float) fom;
	double temp(atanh(fom));
	fdata[inc + 5]          = (float) temp*tab.Cos(pcalc);          // HLA
	fdata[inc + 6]          = (float) temp*tab.Sin_charged(pcalc);  // HLB
	fdata[inc + 7]          = (float) ZERO;                         // HLC
	fdata[inc + 8]          = (float) ZERO;                         // HLD	
	fdata[inc + 9]          = (float) ZERO;  // difference map
	fdata[inc + 10]         = (float) ZERO;  // difference map

	if (phasecomb)
	{
	  // fcomb and phicomb ***NSP
	  fdata[inc + 11]       = (float) (TWO*fom*fobs - dluz*fcalc);
	  fdata[inc + 12]       = (float) (xtal.sf[0].pmodel[r][0]);
	}
      }
      else if (fill)
      {
	fdata[inc + 2]          = (float) (dluz*fcalc);
	fdata[inc + 3]          = (float) (xtal.sf[0].pmodel[r][0]);
	fdata[inc + 4]          = (float) ZERO;
	fdata[inc + 5]          = (float) ZERO;                         // HLA
	fdata[inc + 6]          = (float) ZERO;                         // HLB
	fdata[inc + 7]          = (float) ZERO;                         // HLC
	fdata[inc + 8]          = (float) ZERO;                         // HLD
	fdata[inc + 9]          = (float) ZERO;  // difference map
	fdata[inc + 10]         = (float) ZERO;  // difference map

	if (phasecomb)
	{
	  // fcomb and phicomb ***NSP
	  fdata[inc + 11]       = (float) (dluz*fcalc);
	  fdata[inc + 12]       = (float) (xtal.sf[0].pmodel[r][0]);
	}
      }
    }
    if (outputmtz)
      CMtz::ccp4_lwrefl(MTZOUT, &fdata[0], &colout[0], fdata.size(), r+1);
  }  

  if (outputmtz)
  {
	if (target != "PSAD") {
  	  double overallpluz2(ZERO);
  	  for (unsigned p             = 0; p < xtal.pluz[5].size(); p++)
    	overallpluz2             += xtal.pluz[5][p];
  	  overallpluz2               /= (double)xtal.pluz[5].size();
  	  printf("Overall Luzzati D2 is %.3f\n", overallpluz2);
	}
    
    if (nstop)
    {
      double dnstop((double)nstop);
      string result;
      sprintf(&result[0],"The overall FOM is %.3f\n",
 	      stopfom/dnstop);
      Bp3Result(result.c_str());
    }
    
    CMtz::MtzPut(MTZOUT, " ");
    if (MTZOUT)
      CMtz::MtzFree(MTZOUT);
    if (MTZIN)
      CMtz::MtzFree(MTZIN);
    
  }
  
  if (xtal.verbose                > 1)
    printf("LIKELIHOOD = %f\n",likelihood);
  
  if ( (flagged && check) || (xtal.verbose > 1) )
    printf("Number of reflections with low likelihood values: %u\n", flagged);

  delete[] F;
  delete[] ph;
  
  return likelihood;
}
  
