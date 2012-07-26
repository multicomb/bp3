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

#include "misc.h"
#include "bp3input.h"
#include "crystal.h"
#include "bp3likelihood.h"
#include "minimizer.h"
#include "model.h"
#include <string>
#include "mytimer.h"
Timer Tmain("Main");
Timer Trefine("Refine");
Timer Tmrefine("Minimizer::Refine");
Timer Tlineminimizer("lineminimizer");
Timer Tmore("More");
Timer Tgrad1dim("Gradient1dim");
Timer Tgrad("Grad");
TimerT Tsad("SAD-OUTPUT");
TimerT Tsad1("SAD-NOOUTPUT");
TimerT Tsad2("SAD");
TimerT TinverseGold("InverseGold");
TimerT Tinverse("Inverse");
unsigned long long filter_cnt = 0, filter_cnt_all = 0;

// local functions
void Banner(int, char**);
void GetInput(Model &, Crystal &, Bp3likelihood &, Minimizer &);
void Refine(Bp3likelihood &, Minimizer &);
void Phase(Bp3likelihood &, Minimizer &);
void Check(Bp3likelihood &, Minimizer &);
void Script(Bp3likelihood &);
void Diff(Bp3likelihood &);

int main(int argc, char** argv)
{

  Tmain.start();

  std::string tag;
  tag = Tmain.start("Banner");
  Banner(argc, argv); 
  Tmain.stop(tag);

  // initialize classes for our parameters and data
  Model mdl;
  Crystal xtal;
  Bp3likelihood likelihood(mdl, xtal);
  likelihood.Setcommand(argv[0]);
  Minimizer minimizer(likelihood);
  
  GetInput(mdl, xtal, likelihood, minimizer);

  if (likelihood.Getmode()      == "REFINE")
  {
    tag = Tmain.start("Refine");
    Refine(likelihood, minimizer);
    Tmain.stop(tag);
  }
  else if (likelihood.Getmode() == "DIFF")
  {
    tag = Tmain.start("Diff");
    Diff(likelihood);
    Tmain.stop(tag);
  }
  else if (likelihood.Getmode() == "PHASE")
  {
    tag = Tmain.start("Phase");
    Phase(likelihood, minimizer);
    Tmain.stop(tag);
  }
  else if (likelihood.Getmode() == "CHECK")
  {
    tag = Tmain.start("Check");
    Check(likelihood, minimizer);
    Tmain.stop(tag);
  }

  Tmain.stop();

  Tmain.dump();
  Trefine.dump();
  Tmrefine.dump();
  Tlineminimizer.dump();
  Tmore.dump();
  Tgrad1dim.dump();
  Tgrad.dump();
  Tsad.dump();
  Tsad1.dump();
  Tsad2.dump();
  TinverseGold.dump();
  Tinverse.dump();
  fprintf(stderr, " filter= %llu  all= %llu  frac= %g\n",
      filter_cnt, filter_cnt_all, (double)filter_cnt/(double)filter_cnt_all);
  
  Termination();
}

void Banner(int argc, char** argv)
{
  CCP4Banner(argc,argv);
  printf("                           BP3\n");
  printf(" Multivariate likelihood substructure refinement and phasing\n");
  printf("        http://www.bfsc.leidenuniv.nl/software/bp3\n\n");
}

void GetInput(Model &mdl, Crystal &xtal,
	      Bp3likelihood &like, Minimizer &minimizer)
{
  // gets input and stores it in data structures
  Bp3input input(mdl, xtal, like, minimizer);
}

void Refine(Bp3likelihood &likelihood, Minimizer &minimizer)
{
  Tmrefine.start();
  Trefine.start();
  Trefine.start("Getrefineall");
  if (!likelihood.Getrefineall())
  {
    printf("Refinement cycle\n\n");
    //    likelihood.Setminocc(0.0075);
    //    likelihood.Setmaxluzz(0.5);
    if (likelihood.Gettarget() == "MSRS")
    {
      Trefine.start("Getrefineall::1");
      likelihood.Setrefineocc(false).Setrefineluzzati(true);
      Trefine.stop("Getrefineall::1");
    }
    else
    {
      Trefine.start("Getrefineall::2");
      likelihood.Setrefineocc(true).Setrefineluzzati(false);
      Trefine.stop("Getrefineall::2");
    }
    Trefine.start("Getrefineall::llsetup");
    likelihood.setup(true);
    Trefine.stop("Getrefineall::llsetup");

//    minimizer.Setiterations(std::max((unsigned)10,likelihood.Getnpars()));

    Trefine.start("Getrefineall::refine");
    minimizer.refine();
    Trefine.stop("Getrefineall::refine");

    Trefine.start("Getrefineall::llsetup1");
    likelihood.Setrefineluzzati(true).Setrefineocc(true);
    likelihood.setup();
    Trefine.stop("Getrefineall::llsetup1");

//    minimizer.Setiterations(likelihood.Getnpars());
    Trefine.start("Getrefineall::refine1");
    minimizer.refine();
    Trefine.stop("Getrefineall::refine1");
  }
  Trefine.stop("Getrefineall");

  Trefine.start("setup");
  likelihood.Setrefineluzzati(true).Setrefineocc(true);
  likelihood.Setrefinebfac(true).Setrefinexyz(true); 
  likelihood.setup();
  Trefine.stop("setup");

 // minimizer.Setiterations(3.0*likelihood.Getnpars());
  Trefine.start("refine");
  minimizer.refine();
  Trefine.stop("refine");

  printf("Generating phases and phase probability statistics\n\n");
  Trefine.start("setup1");
  likelihood.setup(false,false,true);
  likelihood.grad(true,true);
  likelihood.Setstats(true);
  likelihood.print();
  Trefine.stop("setup1");

  Trefine.start("script");
  Script(likelihood);
  Trefine.stop("script");

  Trefine.start("Gethand");
  if (likelihood.Gethand())
  {
    // printf("Refinement cycle for other hand (if needed)\n\n");
    likelihood.otherhand();
    likelihood.Setstats(false);

    likelihood.setup(false,true,true);
    /*
    if ( (likelihood.norm()     > 10000.0) && !likelihood.Getmad() )
    {
      likelihood.setup(true);
      minimizer.refine();
      likelihood.setup(false,true,true);
      minimizer.refine();
    }
    else if (likelihood.norm() > 250.0)
      minimizer.refine();
    */

    printf("Generating phases and phase probability statistics for the other hand\n\n");
    likelihood.setup(false,false,true);
    likelihood.grad(true,true);
    likelihood.Setstats(true);
    likelihood.print();

    Script(likelihood);
  }
  Trefine.stop("Gethand");
  Trefine.stop();
  Tmrefine.stop();
}

void Diff(Bp3likelihood &likelihood)
{
  printf("Perform finite difference tests\n\n");
  likelihood.setup(true,true);
  likelihood.finitedifftest();
}

void Phase(Bp3likelihood &likelihood, Minimizer &minimizer)
{
  if (!minimizer.Getsetit())
    minimizer.Setiterations(likelihood.Setiterations());

  if (minimizer.Getiterations())
  {
    if (!likelihood.Getrefineall())
    {
      likelihood.Setminocc(0.025);
      likelihood.Setrefineocc(true).Setrefineluzzati(false);
      likelihood.setup(true);
      minimizer.refine();
      likelihood.Setrefineluzzati(true);
      likelihood.setup();
      minimizer.refine(); 
    }
    likelihood.Setrefinexyz(true);
    likelihood.Setrefinebfac(true);
    likelihood.setup();
    // likelihood.checkatoms();
    //    likelihood.Setsheldrick(false).setup(false,true,true);
    minimizer.refine();
  }
 
  printf("Generating phases and phase probability statistics\n\n");
  likelihood.setup(false,false);
  likelihood.grad(true,true);
  likelihood.Setstats(true);
  likelihood.print();

  Script(likelihood);

  if (likelihood.Gethand())
  {
    //   printf("Refinement cycle for other hand\n\n");
    likelihood.otherhand();
    likelihood.Setstats(false);
    likelihood.setup(false,true,true);
    /*
    if (likelihood.norm() > 250.0)
      minimizer.refine();
    */

    printf("Generating phases and phase probability statistics for the other hand\n\n");
    likelihood.setup(true,false);
    likelihood.grad(true,true);
    likelihood.Setstats(true);
    likelihood.print();

    Script(likelihood);
  }
}

void Check(Bp3likelihood &likelihood, Minimizer &minimizer)
{
  if (!minimizer.Getsetit())
    minimizer.Setiterations(likelihood.Setiterations());
  
  if (minimizer.Getiterations())
  {
    likelihood.Setminocc(0.025);
    likelihood.Setrefineocc(true).Setrefineluzzati(false);
    likelihood.setup(true);
    minimizer.refine();
    likelihood.Setrefineluzzati(true);
    likelihood.setup();
    minimizer.refine();
  }
 
  printf("Generating phases and phase probability statistics\n\n");
  likelihood.setup(false,false);
  likelihood.grad(true,true);
  likelihood.Setstats(true);
  likelihood.print();

  Script(likelihood);

}

void Script(Bp3likelihood &likelihood)
{
  printf("Writing out script for further refinement, crank xml and pdb file\n\n");
  likelihood.writescript();
  likelihood.writexml();
  likelihood.writepdb();
}
