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
#include "minimizer.h"
#include "../mytimer.h"
extern Timer Tmrefine;
extern Timer Tlineminimizer;
extern Timer Tmore;

void Minimizer::refine(const bool savehessian)
{
  std::string tag;
  // Compute function and gradient
  tag = Tmrefine.start("Setpars");
  pars                  = likelihood.Setpars();
  Tmrefine.stop(tag);

  tag = Tmrefine.start("gradient");
  firstfunc             = likelihood.gradient(pars,grad);
  Tmrefine.stop(tag);

  if (grad.size()      <= 0)
  {
    Bp3Warning("Minimizer::refine","no parameters to refine");
    return;
  }

  tag = Tmrefine.start("resize");
  resize();
  func                  = firstfunc; 
  reinit                = 0;
  norm                  = ZERO;
  // initialise Hessian to identity matrix
  negativegamma         = 0;
  Tmrefine.stop(tag);

  if (savehessian)
  {
    tag = Tmrefine.start("savehessian");
    for (unsigned j     = 0; j < grad.size(); j++)
    {
      direction[j]      = ZERO;
      norm             += grad[j]*grad[j];
      for (unsigned k   = 0; k < grad.size(); k++)
        direction[j]   -= hessian[j][k]*grad[k];
    }
    Tmrefine.stop(tag);
  }
  else
  {
    tag = Tmrefine.start("not_savehessian");
    for (unsigned   i   = 0; i < grad.size(); i++)
    {
      for (unsigned j   = 0; j < grad.size(); j++)
        hessian[i][j]   = ZERO;

      hessian[i][i]     = ONE;
      direction[i]      = -grad[i];
      norm             += grad[i]*grad[i];
    }	
    Tmrefine.stop(tag);
  }
  norm                  = sqrt(norm);

  minit                 = std::max(unsigned(grad.size()), (unsigned)5);
  unsigned cyc(maxit);

  if (!setit)
    cyc                 = std::max((unsigned)(2.5*grad.size()), minit);

  printf("Minimum number of iterations: %u\n", minit);
  printf("Maximum number of iterations: %u\n", cyc);

  // Main loop   
  tag = Tmrefine.start("mainloop");
  fprintf(stderr, " >>> cycle= %d \n", cyc);
  for (iter             = 0; iter < cyc ; iter++)
  {
    std::string tag1;
    print();

    tag1 = Tmrefine.start("mainloop::lineminimizer");
    samegrad            = steepest = false;
    oldgrad             = grad;
    oldpars             = pars;
    double oldfunc      = func;
    lineminimizer();
    Tmrefine.stop(tag1);

    tag1 = Tmrefine.start("mainloop::lineminimizer");
    func                = (samegrad) ? newfunc : likelihood.gradient(pars, grad);
    Tmrefine.stop(tag1);

    if ( (reinit        > 5)  ) // || ( ( (oldfunc - func) < DSMALL) && ( (iter < 2) || (iter > (unsigned)2.25*grad.size()) ) ) )   
    {
      Bp3Warning("Minimizer::refine", "minimization without convergence");
      break;
    }


    // test for convergence
    tag1 = Tmrefine.start("mainloop::convtest");
    double partest(ZERO);
    norm                = ZERO;
    for (unsigned j     = 0; j < grad.size(); j++)
    {
      norm             += grad[j]*grad[j];
      double temp       = fabs(pars[j] - oldpars[j])/std::max(ONE, oldpars[j]);
      partest           = (temp > partest) ? temp : partest;
    }
    norm                = sqrt(norm);
    Tmrefine.stop(tag1);

    /*
       if (norm            < ONE)
       break;

       if (norm            < FIVE) 
       break;
       */

    if ( (norm          < normtol) && (iter > minit ) )
      break;

    /* 
       if ( ( (norm        < 50.0)     || (partest < partol) ) && (iter > minit) &&
       (norm          < 7500.0)   && (fabs(oldfunc - func) < 0.01) )
       break;
       */

    // update step direction...
    if (search         == "DFP")
    {
      tag1 = Tmrefine.start("mainloop::dfp");
      dfp();
      Tmrefine.stop(tag1);
    }
    else 
    {
      tag1 = Tmrefine.start("mainloop::bfgs");
      bfgs();
      Tmrefine.stop(tag1);
    }
  }
  Tmrefine.stop(tag);

  if (iter             == cyc)
    Bp3Warning("Minimizer::refine", "maximum number of cycles with no convergence");
  else
    iter++;

  print();
}

void Minimizer::resize()
{
  // Allocate memory for variables needed by the local minimizers
  direction.resize(grad.size());

  hessian.resize(grad.size());
  for (unsigned i = 0; i < grad.size(); i++)
    hessian[i].resize(grad.size());
}

void Minimizer::print() const
{
  //  Print out information detailing the minimization progress 

  printf("Cycle # %u\n", iter+1);
  printf("Log-likelihood      = %10.3f\n", func);
  printf("Log-likelihood gain = %10.3f\n", firstfunc - func);
  printf("Norm of gradient    = %10.3f\n\n", norm);
  fflush(stdout);

  likelihood.print(iter+1);
}

void Minimizer::bfgs()
{
  // updates the bfgs step direction
  vector<double> diffgrad(grad.size(), ZERO);

  for (unsigned j       = 0; j < grad.size(); j++)
  {
    diffgrad[j]         = grad[j] - oldgrad[j];
    direction[j]        = pars[j] - oldpars[j];
  }

  vector<double> hdg(grad.size(), ZERO);

  for (unsigned k       = 0; k < grad.size(); k++)
    for (unsigned j     = 0; j < grad.size(); j++)
      hdg[k]           += hessian[k][j]*diffgrad[j];

  double gamma(ZERO), delta(ZERO);

  for (unsigned j       = 0; j < grad.size(); j++)
  {
    gamma              += diffgrad[j]*direction[j];
    delta              += diffgrad[j]*hdg[j];
  }

  if (gamma             > EPSILON)
  {
    negativegamma       = 0;
    for (unsigned j     = 0; j < grad.size(); j++)
      for (unsigned k   = j; k < grad.size(); k++)
      {
        hessian[j][k]  += ((direction[j]*direction[k]*(delta/gamma + ONE))/gamma - 
            (hdg[j]*direction[k] + direction[j]*hdg[k])/gamma);
        hessian[k][j]   = hessian[j][k];
      }
    for (unsigned j     = 0; j < grad.size(); j++)
    {
      direction[j]      = ZERO;
      for (unsigned k   = 0; k < grad.size(); k++)
        direction[j]   -= hessian[j][k]*grad[k];
    }
  }
  else
  {
    for (unsigned j     = 0; j < grad.size(); j++)
      direction[j]      = -grad[j];
    if (likelihood.Getverbose() > 1)
      Bp3Warning("Minimizer::bfgs", "gamma <= EPSILON");
    negativegamma++;
  }
}

void Minimizer::dfp()
{
  // updates the dfp step direction
  vector<double> diffgrad(grad.size());

  for (unsigned j      = 0; j < grad.size(); j++)
  {
    diffgrad[j]        = grad[j] - oldgrad[j];
    direction[j]       = pars[j] - oldpars[j];
  }

  vector<double> hdg(grad.size());

  for (unsigned k      = 0; k < grad.size(); k++)
  {
    hdg[k]             = ZERO;
    for (unsigned j    = 0; j < grad.size(); j++)
      hdg[k]          += hessian[k][j]*diffgrad[j];
  }

  double gamma(ZERO), delta(ZERO);

  for (unsigned j      = 0; j < grad.size(); j++)
  {
    gamma             += direction[j]*diffgrad[j];
    delta             += diffgrad[j]*hdg[j];
  }

  if ( (gamma          > EPSILON) && (fabs(delta) > EPSILON) )
  {
    for (unsigned j    = 0; j < grad.size(); j++)
      for (unsigned k  = j; k < grad.size(); k++)
      {
        hessian[j][k] += (direction[j]*direction[k]/gamma -
            hdg[j]*hdg[k]/delta);
        hessian[k][j]  = hessian[j][k];
      }
  }
  else
    if (likelihood.Getverbose() > 1)
      Bp3Warning("Minimizer::dfp", "gamma <= EPSILON or fabs(delta) <= EPSILON");

  for (unsigned j      = 0; j < grad.size(); j++)
  {
    direction[j]       = ZERO;
    for (unsigned k    = 0; k < grad.size(); k++)
      direction[j]    -= hessian[j][k]*grad[k];
  }
}

void Minimizer::lineminimizer()
{
  Tlineminimizer.start("Main");
  lineconverge        = steepest = false;
  nfe                 = 0;

  std::string tag;

  tag = Tlineminimizer.start("Main::dot_prod");
  stepmax             = ONE;   // ***NSP
  dfunc               = dot_prod(grad,direction);
  Tlineminimizer.stop(tag);


  if (dfunc          >= ZERO)
  {
    if (likelihood.Getverbose() > 1)
      Bp3Warning("Minimizer::lineminimizer", 
          "switching to steepest descent");
    for (unsigned i   = 0; i < grad.size(); i++)
      direction[i]    = -grad[i];
  }  

  if (likelihood.Getverbose())
  {
    printf("Line minimization\n");
    printf("Initial step size:      %g\n", step);
  }

  if (line           == "MORE")
  {
    tag = Tlineminimizer.start("Main::MORE");
    if (!more())
    {
      if ( fbest      < func ) 
      {
        step          = stepbest;
        newfunc       = fbest;
        dnewfunc      = dfbest;
      }
      else
      {
        newfunc       = func;
        dnewfunc      = dfunc;
        step          = ZERO;
      }
    }
    else 
      lineconverge    = true;
    if (fabs(newfunc  - likelihood.Getfunction()) < DSMALL)
      samegrad        = true;
    Tlineminimizer.stop(tag);
  }
  else
  {
    tag = Tlineminimizer.start("Main::notMORE");
    newfunc           = likelihood.function1dim(step,pars,direction); nfe++;

    if (newfunc       > (func + alpha*step*dfunc))
    { 
      double diff     = fabs((newfunc - func)/func);
      step           /= std::max(ONE,std::min(diff*10.0, 1000.0));
      newfunc         = brent(step);
      lineconverge    = (step > DSMALL);
    }
    Tlineminimizer.stop(tag);
  }

  if (step            < SMALLESTD)  
  {
    tag = Tlineminimizer.start("Main::step_small");
    steepest          = true;
    reinit++;
    if (reinit        > 1)
      reinitializehessian();
    likelihood.function1dim(0.0,pars,direction);
    Tlineminimizer.stop(tag);
  }
  else
  {
    tag = Tlineminimizer.start("Main::step_else");
    if (fabs(newfunc  - likelihood.Getfunction()) > 0.0001)
      likelihood.function1dim(step,pars,direction);
    pars              = likelihood.Setpars();
    //    reinit            = 0;
    Tlineminimizer.stop(tag);
  }

  if (likelihood.Getverbose())
  {
    printf("Final step size:        %g\n", step);
    printf("Function evaluations:   %u\n\n", nfe);
  }

  tag = Tlineminimizer.start("Main::check_shift");
  checkshift();
  Tlineminimizer.stop(tag);

  if (step            < EPSILON)
    step              = 0.00001;
  Tlineminimizer.stop("Main");
}

bool Minimizer::more()
{
  Tmore.start("Main");
  // Reference:
  // J. J. More' and D. J. Thuente,
  //  Line search algorithms with guaranteed sufficient decrease,
  //  ACM Transactions on  Mathematical Software, 20, (1994), 286--307

  bool interval(false), converge(false);
  double width(stepmax-stepmin), width1(width*TWO);
  double gtest(alpha*dfunc), xtol(0.001);
  double stmin(ZERO);
  double stmax(ONE);
  double fx(func), fy(func), gx(dfunc), gy(dfunc);
  double stx(ZERO), sty(ZERO);
  static vector<double> flast;
  flast.clear();
  fbest           = func;
  dfbest          = dfunc;
  stepbest        = ZERO;
  bracket         = false;

  for (unsigned i = 0; i < maxlineit; i++)
  {
    std::string tag;

    tag = Tmore.start("Main::gradient1dim");
    newfunc       = likelihood.gradient1dim(dnewfunc, grad, step, pars, direction);
    Tmore.stop(tag);
    flast.push_back(newfunc);
    nfe++;

    tag = Tmore.start("Main::misc");
    if (newfunc   < fbest)
    {
      stepbest    = step;
      fbest       = newfunc;
      dfbest      = dnewfunc;
    }

    double ftest  = func + step * gtest;
    if (!interval && newfunc <= ftest && dnewfunc >= ZERO)
      interval    = true;

    if (bracket && (step <= stmin || step >= stmax))
    {
      if (likelihood.Getverbose() > 1)
        Bp3Warning("Minimizer::more", "Roundings errors prevent minimization");
      steepest    = true;
      break;
    }
    if (step     == stepmax && newfunc <= ftest && dnewfunc <= gtest)
    {
      converge    = true;
      break;
    }    
    if (bracket && stmax - stmin <= xtol * stmax)
    {
      converge    = true;
      break;
    }        
    //  Test for convergence.
    if (newfunc  <= ftest && fabs(dnewfunc) <= beta * (-dfunc))
    {
      converge    = true;
      break;
    }

    if (i         > 2)
      if ( (fabs(flast[i] - flast[i-1]) < 1.0e-6) &&
          (fabs(flast[i] - flast[i-2]) < 1.0e-6) &&
          (fabs(flast[i] - flast[i-3]) < 1.0e-6)   )
      {
        Bp3Warning("Minimizer::more", "Function not decreasing anymore");
        break;
      }
    Tmore.stop(tag);

    if (!interval && newfunc <= fx && newfunc > ftest) 
    {
      tag = Tmore.start("Main::choosestep");
      double fm   = newfunc - step * gtest;
      double fxm  = fx - stx * gtest;
      double fym  = fy - sty * gtest;
      double gm   = dnewfunc - gtest;
      double gxm  = gx - gtest;
      double gym  = gy - gtest;
      choosestep(stx, fxm, gxm, sty, fym, gym, fm, gm, stmin, stmax);
      fx          = fxm + stx * gtest;
      fy          = fym + sty * gtest;
      gx          = gxm + gtest;
      gy          = gym + gtest;
      Tmore.stop(tag);
    } 
    else
    {
      tag = Tmore.start("Main::choosestep1");
      choosestep(stx, fx, gx, sty, fy, gy, newfunc, dnewfunc, stmin, stmax);
      Tmore.stop(tag);
    }

    tag = Tmore.start("Main::misc2");
    if (bracket) 
    {
      if (fabs(sty - stx) >= width1 * TWOTHIRD)
        step      = stx + (sty - stx) * HALF;
      width1      = width;
      width       = fabs(sty - stx);
    }
    if (bracket)  
    {
      stmin       = std::min(stx,sty);
      stmax       = std::max(stx,sty);
    } 
    else 
    {
      stmin       = step + (step - stx) * 1.1;
      stmax       = step + (step - stx) * FOUR;
    }
    step          = std::min(std::max(step,stepmin), stepmax);
    if ((bracket && (step <= stmin || step >= stmax)) || 
        (bracket && stmax - stmin <= xtol * stmax))
      step        = stx;
    Tmore.stop(tag);
  }
  Tmore.stop("Main");
  return converge;
} 

void Minimizer::choosestep(double &stx,  double &fx, double &dx, 
    double &sty,  double &fy, double &dy,
    double &fp, double &dp,
    const double stmin, const double stmax)
{

  double stepf(ZERO), stepc(ZERO), stepq(ZERO);

  double sgnd        = dp * (dx / fabs(dx));

  if (fp             > fx)
  {
    double theta     = (fx - fp) * THREE / (step - stx) + dx + dp;
    double s         = std::max(std::max(fabs(theta),fabs(dx)), fabs(dp));
    double gamma     = s * sqrt((theta*theta)/(s * s) - dx / s * (dp / s));
    if (step         < stx) 
      gamma         *= -ONE;
    double p         = gamma - dx + theta;
    double q         = gamma - dx + gamma + dp;
    double r         = p / q; 
    stepc            =  stx + r * (step - stx);
    stepq            = stx + dx / ((fx - fp) / (step - stx) + dx) / TWO * (step - stx);
    stepf            = ( (fabs(stepc - stx) < fabs(stepq - stx)) ? 
        stepc : stepc + (stepq - stepc)/TWO);
    bracket          = true;
  } 
  else if (sgnd      < ZERO) 
  {
    double theta     = (fx - fp) * THREE / (step - stx) + dx + dp;
    double s         = std::max(std::max(fabs(theta), fabs(dx)), fabs(dp));
    double gamma     = s * sqrt(theta*theta/(s*s) - dx / s * (dp / s));
    if (step         > stx) 
      gamma         *= -ONE;
    double p         = gamma - dp + theta;
    double q         = gamma - dp + gamma + dx;
    double r         = p / q;
    stepc            = step + r * (stx - step);
    stepq            = step + dp / (dp - dx) * (stx - step);
    stepf            = (fabs(stepc - step) > fabs(stepq - step)) ? stepc : stepq;
    bracket          = true;
  } 
  else if (fabs(dp)  < fabs(dx)) 
  {
    double theta     = (fx - fp) * THREE / (step - stx) + dx + dp;
    double s         = std::max(std::max(fabs(theta),fabs(dx)) ,fabs(dp));
    double gamma     = s * sqrt((std::max(ZERO,
            theta*theta/(s*s) - dx / s * (dp / s))));
    if (step         > stx)
      gamma          = -ONE;
    double p         = gamma - dp + theta;
    double q         = gamma + (dx - dp) + gamma;
    double r         = p / q;
    if (r            < ZERO && gamma != ZERO)  
      stepc          = step + r * (stx - step);
    else if (step    > stx) 
      stepc          = stmax;
    else
      stepc          = stmin;
    stepq            = step + dp / (dp - dx) * (stx - step);
    if (bracket) 
    {
      stepf          = (fabs(stepc - step) < fabs(stepq - step)) ? stepc : stepq;
      if (step       > stx) 
        stepf        = std::min(step + (sty - step) * TWOTHIRD,stepf);
      else 		 	
        stepf        = std::max(step + (sty - step) * TWOTHIRD,stepf);
    } 
    else
    {
      stepf          = (fabs(stepc - step) > fabs(stepq - step) ) ? stepc : stepq;
      stepf          = std::max(std::min(stmax,stepf), stmin);
    }
  } 
  else 
  {
    if (bracket) 
    {
      double theta   = (fp - fy) * THREE / (sty - step) + dy + dp;
      double s       = std::max(std::max(fabs(theta),fabs(dy)),fabs(dp));
      double gamma   = s * sqrt(theta*theta/(s*s) - dy / s * (dp / s));
      if (step       > sty) 
        gamma        = -ONE;
      double p       = gamma - dp + theta;
      double q       = gamma - dp + gamma + dy;
      double r       = p / q;
      double stepc   = step + r * (sty - step);
      stepf          = stepc;
    } 
    else if (step    > stx) 
      stepf          = stmax;
    else
      stepf          = stmin;    
  }
  if (fp             > fx) 
  { 
    sty              = step;
    fy               = fp;
    dy               = dp;
  } 
  else 
  {
    if (sgnd         < ZERO) 
    {
      sty            = stx;
      fy             = fx;
      dy             = dx;
    }
    stx              = step;
    fx               = fp;
    dx               = dp;
  }
  step               = stepf;
}

double Minimizer::brent(double &x)
{
  double a          = stepmin;
  double b          = stepmax;
  double v(ZERO);
  double w          = v;
  x                 = v;
  double fx         = func;
  double fv         = fx;
  double fw         = fx;

  double d(ZERO), e(ZERO), u(ZERO), fu(ZERO);

  for (unsigned i   = 0; i < 100; i++) 
  {
    double xm       = HALF*(a + b);
    double eps      = 0.2*fabs(x) + 1.0e-10;
    double tol      = eps*TWO;

    if (fabs(x-xm) <= tol - (b - a)*HALF) 
      break;

    bool change(true);

    if (fabs(e)     > eps) 
    {
      double r      = (x - w)*(fx - fv);
      double q      = (x - v)*(fx - fw);
      double p      = (x - v)*q - (x - w)*r;
      q             = (q - r)*TWO;
      if (q         > ZERO)
        p          *= -ONE;
      q             = fabs(q);
      double temp   = e;
      e             = d;
      change        = (fabs(p) >= fabs(q*HALF*temp) ||
          p <= q * (a - x) || p >= q * (b - x));
      if (!change)
      {
        d           = p / q;
        u           = x + d;
        if (u - a   < tol || b - u < tol) 
          d         = dsign(eps, xm - x);
      }
    }

    if (change)
    {
      if (x        >= xm) 
        e           = a - x;
      else 
        e           = b - x;
      d             = e * 0.381966;
    }

    if (fabs(d)    >= eps) 
      u             = x + d;
    else 
      u             = x + dsign(eps, d);
    fu              = likelihood.function1dim(u,pars,direction); nfe++;
    if (fu         <= fx) 
    {
      if (u        >= x) 
        a           = x;
      else 
        b           = x;
      v             = w;
      fv            = fw;
      w             = x;
      fw            = fx;
      x             = u;
      fx            = fu;
    } 
    else 
    {
      if (u         < x) 
        a           = u;
      else 
        b           = u;
      if (fu       <= fw || w == x) 
      {
        v           = w;
        fv          = fw;
        w           = u;
        fw          = fu;
      } 
      else if (fu  <= fv || v == x || v == w) 
      {
        v           = u;
        fv          = fu;
      }
    }
  }
  return fx;
}

void Minimizer::checkshift()
{
  // check to see if the shift were damped and if they are
  // modify the direction accordingly

  unsigned changed(0);

  for (unsigned i    = 0; i < pars.size(); i++)
    if (fabs((pars[i] - oldpars[i]) - direction[i]*step) > 0.00001)
      if (step       > SMALLESTD)
      {
        direction[i] = (oldpars[i] - pars[i])/step;
        changed++;
      }

  if (changed)
    printf("direction changed for %u parameter(s)\n",changed);
}
