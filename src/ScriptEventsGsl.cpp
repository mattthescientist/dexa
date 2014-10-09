// ScriptEventsGsl.cpp
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
#include "ScriptReader.h"
#include "MgstFcnGsl.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


//------------------------------------------------------------------------------
// Command : startfit <x> <y> <z>
// Collates all the script information (which must have been given before this
// command), constructs the appropriate Minuit fitting objects, and attempts to
// optimise the user's fit parameters based on their initial guesses. Once the
// fit is complete, the parameter values that Minuit considers to be optimal are
// transferred to the spectra stored in Spectra, transforming them into the
// 'fitted' spectra. Since the Minuit objects are defined within the scope of
// this function, no further analysis of the fit may be performed once this
// function terminates.
//
int startfit (stringstream &args, vector<struct spectrum> &Spectra) {

  int Iterations, Verbose;
  MgstFcn theFcn (&Spectra);
  if (args >> Iterations >> Verbose) {
    theFcn.iterations (Iterations);
  } else { 
    args >> Verbose; 
    Iterations = DEF_MAX_ITERATIONS;
  }

  cout << 
    "========================================================================="
    << endl << "Starting DiffXAS fit with GSL ..." << endl << endl;
  MgstFcn GslFit (&Spectra);

  // Prepare the parameter and error arrays for GSL.
  vector<double> UserParams;
  vector<double> ParamErrors;
  vector<double> MinusLimits;
  vector<double> PlusLimits;
  vector<string> ParamNames;
  if (Verbose >= MAX_VERBOSE) {
	  cout << "...creating parameter arrays" << endl;
  }
  createParameterArrays (UserParams, MinusLimits, PlusLimits, ParamErrors, ParamNames, Spectra, Verbose);
  if (Verbose >= MAX_VERBOSE) {
  	  cout << "...running GslFit ()" << endl;
  }
  // Execute the fit
  GslFit (UserParams);
  cout << "Fitting Complete" << endl << endl;

  // Unlike MINUIT, there is no need to apply the optimised parameter values to
  // the spectra since that is done implictly by GSL. Just print the results.
  printParams (cout, Spectra);
  return SCRIPT_NO_ERROR;
}
