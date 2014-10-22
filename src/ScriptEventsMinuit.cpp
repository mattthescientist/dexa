// ScriptEventsMinuit.cpp
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
#include "MgstFcnMinuit.h"
#include "Minuit/FunctionMinimum.h"
#include "Minuit/MnUserParameters.h"
#include "Minuit/MnMigrad.h"
#include "Minuit/MnMinimize.h"
#include "Minuit/MnSimplex.h"
#include "Minuit/MnMachinePrecision.h"
#include "Minuit/MnMinos.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <csignal>

#define MINUIT_STRATEGY 0   // 0 for basic fitting. Higher values fit noise.
#define MINUIT_ERR_DEF  1.0 // 1 for Chi^2 fitting as described in Minuit manual


//------------------------------------------------------------------------------
// printFitStrategy : Given an MnStrategy object, print its various properties
// on the standard output. This is only used for high verbose levels.
//
void printFitStrategy (MnStrategy &Strategy) {
  cout << "Gradient Cycles    : " << Strategy.gradientNCycles() << endl;
  cout << "Gradient Step Tol  : " << Strategy.gradientStepTolerance() << endl;
  cout << "Gradient Tolerance : " << Strategy.gradientTolerance() << endl;
  cout << endl;
  cout << "Hessian Cycles    : " << Strategy.hessianNCycles() << endl;
  cout << "Hessian Step Tol  : " << Strategy.hessianStepTolerance() << endl;
  cout << "Hessian G2 Tol    : " << Strategy.hessianG2Tolerance() << endl;
  cout << "Hessian Grad Cyc  : " << Strategy.hessianGradientNCycles() << endl;
}


//------------------------------------------------------------------------------
// applyOptimisedParameters : Extracts the parameter values from a 
// FunctionMinimum object passed in as arg1, and applies them to their 
// respective spectra, stored in the spectrum vector at arg2.
//
void applyOptimisedParameters (FunctionMinimum &min, 
  vector<struct spectrum> &Spectra) {
  
  for (unsigned int i = 0; i < Spectra.size (); i ++) {
    for (unsigned int j = 0; j < Spectra[i].SubStructure.size (); j ++) {
      for (unsigned int k = 0; k < Spectra[i].SubStructure[j].Crystal.size (); k ++) {
      (Spectra[i].SubStructure[j].Spectrum.*Spectra[i].SubStructure[j].Crystal[k].ParamSet)
        (min.userState().value (Spectra[i].SubStructure[j].Crystal[k].Index));
      }
      for (unsigned int k = 0; k < Spectra[i].SubStructure[j].Path.size (); k ++) {
        (Spectra[i].SubStructure[j].Spectrum.*Spectra[i].SubStructure[j].Path[k].ParamSet)
          (Spectra[i].SubStructure[j].Path[k].PathNum, min.userState().value 
          (Spectra[i].SubStructure[j].Path[k].Index));
      }
    }
  }
}

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
int startfit(stringstream &args, vector<struct spectrum> &Spectra){

  // Load in the command line parameters. Return an error if they're bad.
  int Iterations, Verbose;
  MgstFcn theFcn (&Spectra);
  args >> Verbose;
  if (!args.eof ()) {
    Iterations = Verbose;
    theFcn.iterations (Iterations);
    args >> Verbose;
  } else { Iterations = DEF_MAX_ITERATIONS; }
  if (!args.eof ()) { return SCRIPT_BAD_ARGS; }

  // Inform the user that the fitting is about to start 
  cout << endl;
  if (Verbose >= MIN_VERBOSE) {
    cout << "================================================================================" << endl;
  }
  if (Verbose >= MED_VERBOSE) {
    cout << "Maximum allowed iterations : " << Iterations << endl;
  }
  
//  (void) signal (SIGINT, catchSigInt);
  
  // Prepare the parameter and error arrays for MINUIT.
  vector<double> UserParams;
  vector<double> ParamErrors;
  vector<double> MinusLimits;
  vector<double> PlusLimits;
  vector<string> ParamNames;
  double NumIndepPoints = createParameterArrays (UserParams, MinusLimits, 
    PlusLimits, ParamErrors, ParamNames, Spectra, Verbose);
  MnUserParameters Parameters (UserParams, ParamErrors);
  for (unsigned int i = 0; i < UserParams.size (); i ++) {
    if (MinusLimits[i] != PlusLimits[i]) {
      Parameters.setLimits (i, MinusLimits[i], PlusLimits[i]);
    }
  }

  // Define the fitting strategy and print its details on high verbose levels
  MnStrategy Strategy (MINUIT_STRATEGY);
  if (Verbose >= MAX_VERBOSE) printFitStrategy (Strategy);

  // Perform the fit then apply the optimised parameter values to the spectra
  cout << "Starting DiffXAS fit with Minuit ..." << flush;
  MnMigrad MiGrad (theFcn, Parameters, Strategy);
  theFcn.setErrorDef(MINUIT_ERR_DEF);
  theFcn.setVerbose (Verbose);
  FunctionMinimum min = MiGrad ();
  applyOptimisedParameters (min, Spectra);
  cout << " done" << endl;

  // Now assess the parameter errors using Minuit's MINOS routines. This step
  // isn't possible if Minuit failed to return a valid minimum. In this case
  // just print a warning message to the standard output.
  if (min.isValid ()) {
    cout << "Calculating 1-sigma parameter errors ..." << flush;
    MnMinos minos (theFcn, min);
    ParamErrors.clear ();
    for (unsigned int i = 0; i < UserParams.size (); i ++) {
      ParamErrors.push_back (minos(i).second);
    }
    cout << " done" << endl;
  } else {
    cout << "Warning : Minuit says the function minimum is invalid." << endl;
    cout << "          Unable to calculate parameter errors." << endl;
  }
  
  // Print details of the fit, including the optimised parameter values and
  // errors, to the standard output. On high verbose levels extra information 
  // is given.
  if (Verbose >= MED_VERBOSE) {
    cout << endl;
    cout << "iterations used        = " << min.nfcn () << endl;
    cout << "distance to minimum    = " << min.edm () << endl;
    cout << "num independent points = " << NumIndepPoints << endl;
    cout << "num parameters fitted  = " << UserParams.size () << endl;
  }
  cout << "chi-squared         = " << min.fval () << endl;
  cout << "reduced chi-squared = " 
    << min.fval () / (NumIndepPoints - double(UserParams.size ())) << endl;
  cout << endl << "Optimised parameter values are:" << endl;
  for (unsigned int i = 0; i < UserParams.size (); i ++) {
    cout << "  " << ParamNames[i] << " = " << min.userState().value(i)
      << " +/- " << ParamErrors[i] << endl;
  }
  cout << endl;

  // No errors encountered, so return with SCRIPT_NO_ERROR
  return SCRIPT_NO_ERROR;
}

/*void catchSigInt (int sig) {
  cout << endl << "Terminating fit. Attempting to finish processing script. Press Ctrl-C again to abort." << endl;
  (void) signal (SIGINT, SIG_DFL);
}*/
