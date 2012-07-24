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
#include "bp3input.h"
#include "ccp4_parser.h"                  // Peter Briggs' parser in c libraries
#include "likelihood.h"

Bp3input::Bp3input(Model &model, Crystal &xtl, Bp3likelihood &likelihood,
		   Minimizer &minimizer)
  : Input(model, xtl, likelihood, minimizer), like(likelihood)
{
  
  // parse
  CCP4parse();

  CCP4::ccp4_parse_end(parser); // clean up parser array
  
  check();

  print();
}

void Bp3input::check() 
{
  // checks if input is okay

  Input::check();

  for (unsigned d = 0; d < xtal.sf.size(); d++)
  {
    if (like.Getmode()  == "CHECK")
      if (!xtal.sf[d].Getnbins())
        xtal.sf[d].Setusernbins(xtal.sf[d].GetMINBINS());

    if ( (xtal.sf[d].Getsdatap() == "") || (xtal.sf[d].Getsdevp() == "") )
      Bp3Error("Bp3Input::check", "Dataset declared without F and/or SIGF columns");
    if ( (xtal.sf[d].Gettype()  != "AMPLITUDE") && (like.Getmode() != "OUTPUTHKL") &&
	 (like.Getmode()        != "OUTPUTSCALEPACK") )
      Bp3Error("Bp3likelihood::setup","input amplitudes only");    
  }
 
  if ( !((like.Getmode() == "INTENSITY") || (like.Getmode() == "OUTPUTHKL") || 
	(like.Getmode()  == "OUTPUTSCALEPACK") ) )
  {  
    if (mdl.site.size() <= 0)
    {
      Bp3Warning("Bp3Input::check", "no atoms have been defined");
      Bp3Warning("Bp3Input::check", "BP3 will output intensities");
      if (like.Getmode() != "OUTPUTSCALEPACK") 
        like.Setmode("INTENSITY");
    }

    if ( (xtal.sf.size()    == 1) && (xtal.sf[0].Getsdatam() == "") )
    {
      Bp3Warning("Bp3Input::check", "only one nonanomalous data set is given");
      Bp3Warning("Bp3Input::check", "BP3 will output intensities"); 
      if (like.Getmode() != "OUTPUTSCALEPACK") 
        like.Setmode("INTENSITY");
    }
  }
}

void Bp3input::CCP4parse()
{
  // The setup of reading keyword input and it follows pretty much
  // exactly Peter Briggs' documentation of his parser.

  // read lines from stdin until END or EOF is reached
  bool cont       = true;
  char line[301];

  while (cont)
  {
    // Blank line before calling ccp4_parser implies reading from stdin
    line[0]       = '\0';
    unsigned ntok = CCP4::ccp4_parser(line, 300, parser, 1);
    
    if (ntok      < 1)                                            // reached EOF
      break;
    else   // perform keyword analysis and interpretation
    {
      if (Input::xtalkeywords())
	continue;
      else if (Input::datakeywords())
	continue;
      else if (Input::modelkeywords())
	continue;
      else if (Input::minimizerkeywords())
	continue;
      else if (Input::likelihoodkeywords())
	continue;
      else if (bp3keywords())
	continue;
      else if (CCP4::ccp4_keymatch("END", parser->keyword))
	break;
      else
	Bp3Error("CCP4parse","Keyword not recognised");
    }
  }  
}

bool Bp3input::bp3keywords()
{
  if (CCP4::ccp4_keymatch("REFI", parser->keyword))       // mode options
    refikeyword();
  else if (CCP4::ccp4_keymatch("PHAS", parser->keyword))            
    phaskeyword();
  else if (CCP4::ccp4_keymatch("CHEC", parser->keyword))            
    checkeyword();
  else if (CCP4::ccp4_keymatch("HAND", parser->keyword))            
    handkeyword();
  else if (CCP4::ccp4_keymatch("DIFF", parser->keyword))            
    diffkeyword();
  else if (CCP4::ccp4_keymatch("SHEL", parser->keyword))            
    shelkeyword();
  else if (CCP4::ccp4_keymatch("FHOU", parser->keyword))            
    fhoukeyword();
  else if (CCP4::ccp4_keymatch("THRE", parser->keyword))            
    threkeyword();  
  else
    return false;
  
  return true;
}

void Bp3input::refikeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 1)
    Bp3Warning("Bp3input::refikeyword", "Extra text in REFIne keyword will be ignored");
  
  like.Setmode("REFINE");
}

void Bp3input::phaskeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Bp3input::phaskeyword", "Extra text in PHASe keyword will be ignored");
  
  like.Setmode("PHASE");
  like.Setmaxocc(HALF).Setmaxluzz(ONE);

  if (parser->ntokens == 2)
    if (parser->token[1].isnumber)
      mini.Setiterations((unsigned)parser->token[1].value);  
    else
      Bp3Error("Bp3input::phaskeyword", "A numerical value is expected after PHAS");
}

void Bp3input::checkeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Bp3input::checkeyword", "Extra text in CHECk keyword will be ignored");
  
  like.Setmode("CHECK");

  if (parser->ntokens == 2)
    if (parser->token[1].isnumber)
      mini.Setiterations((unsigned)parser->token[1].value);  
    else
      Bp3Error("Bp3input::checkeyword", "A numerical value is expected after CHECk");
}

void Bp3input::handkeyword() const
{
  Bp3Warning("Bp3input::handkeyword", "This keyword is no longer needed - it is a default!");

  // check to see if we have no more than one token
  if (parser->ntokens > 1)
    Bp3Warning("Bp3input::handkeyword", "Extra text in HAND keyword will be ignored");
  
  like.Sethand(1);
}

void Bp3input::diffkeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Bp3input::diffkeyword", "Extra text in DIFF keyword will be ignored");

  like.Setmode("DIFF");
  
  if (parser->ntokens == 2)
    if (parser->token[1].isnumber)
      like.Setdifftol((double)parser->token[1].value);
    else
      Bp3Error("Bp3input::diffkeyword", "A numerical value is expected after DIFF");
}

void Bp3input::shelkeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens > 2)
    Bp3Warning("Bp3input::shelkeyword", "Extra text in SHELdrick keyword will be ignored");
  if (parser->ntokens == 2)
    if (parser->token[1].isnumber)
      like.Setshelthres((unsigned)parser->token[1].value);
    else
      Bp3Error("Bp3input::shelkeyword", "A numerical value is expected after SHEL");
  else
    like.Setshelthres(0);
}

void Bp3input::fhoukeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 1)
    Bp3Warning("Bp3input::fhoukeyword", "Extra text in FHOUt keyword will be ignored");

  like.Setoutputhcalc(true);
}
 
void Bp3input::threkeyword() const
{
  // check to see if we have no more than one token
  if (parser->ntokens  > 2)
    Bp3Warning("Bp3input::threkeyword", "Extra text in THRE keyword will be ignored");

  if (parser->token[1].isnumber)
    like.Setthreshold((double)parser->token[1].value);
  else
    Bp3Error("Bp3input::threkeyword", "A numerical value is expected for the THREshold");
}
