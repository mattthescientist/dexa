// MgstFit.cpp
// Copyright (C) 2007-2009 Matthew Ruffoni
//
// This file is part of DEXA.
//
// DEXA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DEXA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with DEXA.  If not, see <http://www.gnu.org/licenses/>.
//
// The code may be compiled using 'make'. Just enter any one of the following:
//
//  > make          Makes both the GSL and Minuit versions of FitChi
//  > make gsl      Makes only the GSL version of FitChi
//  > make minuit   Makes only the Minuit version of FitChi
//  > make clean    Run after calling any of the above to remove sre obj files
//
// Alternatively, in the absence of 'make', the following commands can be used
// For the GSL Scaled Levenberg-Marquart back-end, enter:
//
//   g++ MgstSpectrum.cpp DiffXasSpectrum.cpp MgstTensor.cpp PathFinder.cpp MgstFcnGsl.cpp MgstFitGsl.cpp spectrum.cpp -lgsl -lgslcblas
//     
// And for the Minuit MiGrad back-end, enter:
//
//   g++ MgstSpectrum.cpp DiffXasSpectrum.cpp MgstTensor.cpp PathFinder.cpp MgstFcnMinuit.cpp MgstFitMinuit.cpp spectrum.cpp -lgsl -lgslcblas -llcg_Minuit

#include "ScriptReader.h"
#include <vector>
#include <cstring>

#define MIN_CMD_LINE_PARAMS 2   // i.e. the binary plus one argument
#define MAX_CMD_LINE_PARAMS 3   // i.e. the above plus a verbose switch

#define SCRIPT_ARG  1
#define VERBOSE_ARG 2

#define VERSION_TAG "dexa 1.0.1"

void printHelp () {
  cout << "Syntax:\n\
  dexa <input script> [OPTIONS]...\n\n\
Provides quantitative analysis of DiffXAS spectra.\n\
See \"http://sourceforge.net/projects/dexa\" for more information.\n" << endl;
  cout << "Options:\n\
  <input script> : specifies the name of a DEXA script to be processed. If this\n\
                   file is not present in the current directory, its path must\n\
                   also be given.\n\n\
  -v[0|1|2|3] : Verbose output. By default, the output level is 1, which is\n\
                equivalent to '-v1'. Specifying '-v' without an explicit level\n\
                will raise the output to level 2. Each levels does the following\n\n\
                -v0 : Minimal output. Prints the fit parameters and errors, and\n\
                      saves only the fully-processed spectra.\n\
                -v1 : The same output as '-v0' plus some extra files, typically\n\
                      showing intermediate levels of various calculations.\n\
                -v2 : Provides the full user-level output. Equivalent to '-v'.\n\
                -v3 : Provides all the above plus debugging information.\n\n";
}
  

int main (int argc, char* argv[]) {
  std::vector <struct spectrum> theSpectra;
  int Verbose = MIN_VERBOSE;
  
  cout << VERSION_TAG << endl;
  // The user has run DEXA by simply specifying a script. Turn verbose mode off.
  if (argc == MIN_CMD_LINE_PARAMS) {
    readScript (argv[SCRIPT_ARG], theSpectra, Verbose);
    
  // The user has supplied a verbose argument. Set verbose level accordingly.
  } else if (argc == MAX_CMD_LINE_PARAMS) {
   Verbose = ERR_VERBOSE;
    if (strcmp (argv[VERBOSE_ARG],"-v0") == 0) Verbose = NO_VERBOSE;
    else if (strcmp (argv[VERBOSE_ARG],"-v1") == 0) Verbose = MIN_VERBOSE;
    else if (strcmp (argv[VERBOSE_ARG],"-v")  == 0) Verbose = MED_VERBOSE;
    else if (strcmp (argv[VERBOSE_ARG],"-v2") == 0) Verbose = MED_VERBOSE;
    else if (strcmp (argv[VERBOSE_ARG],"-v3") == 0) Verbose = MAX_VERBOSE;
    if (Verbose != ERR_VERBOSE) {
      readScript (argv[SCRIPT_ARG], theSpectra, Verbose);
    } else {
      cout << "Incorrect command line parameters. ";
      printHelp ();
      return SCRIPT_SYNTAX_ERROR;
    }      

  // No match found so print some help and return a syntax error
  } else {
    cout << "Incorrect command line parameters. ";
    printHelp ();
    return SCRIPT_LOAD_ERROR;
  }
  return SCRIPT_NO_ERROR;
}
