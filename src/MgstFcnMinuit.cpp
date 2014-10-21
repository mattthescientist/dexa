// MgstFcnMinuit.cpp
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
#include "MgstFcnMinuit.h"
#include <string>
#include <iostream>
#include <cmath>
#include <cassert>
#include <sstream>

using namespace::std;

// This class inherits Minuit's FCNBase class and so provides the fitting 
// function (Fcn) for MiGrad. There is just one member function here: the
// overloaded operator (). When called, it will insert the fit parameters into
// the spectrum specified in the class constructor and then calculate the sum of
// the differences squared between this theory spectrum and the experimental
// data. This value is then returned so that MiGrad may then continue with its
// fitting process.
//
// Presently two spectra may be fitted at once hence the two Chi2 loops below.

double MgstFcn::operator() (const vector<double> &par) const {
  
  static int Iteration = 0;
  Iteration ++;
  if (Iteration > MaxIterations) return 0.0;
  
  
  // Use the Index property of each fit parameter in theSpectra to apply the new
  // parameter values to the actual MgstSpectrum objects.
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    for (unsigned int k = 0; k < theSpectra -> at(i).SubStructure.size (); k ++) {

      // Apply the parameters relating to specific scattering paths first...
      for (unsigned int j = 0; j < theSpectra -> at (i).SubStructure[k].Path.size (); j ++) {
        (theSpectra -> at (i).SubStructure[k].Spectrum.*theSpectra -> at (i).SubStructure[k].Path[j].ParamSet)
          (theSpectra -> at (i).SubStructure[k].Path[j].PathNum, par[theSpectra -> at (i).SubStructure[k].Path[j].Index]);
      }

      // ... and then the parameter relating to the crystal as a whole
      for (unsigned int j = 0; j < theSpectra -> at (i).SubStructure[k].Crystal.size (); j ++) {
        (theSpectra -> at (i).SubStructure[k].Spectrum.*theSpectra -> at (i).SubStructure[k].Crystal[j].ParamSet)
          (par[theSpectra -> at (i).SubStructure[k].Crystal[j].Index]);
      }
    }
  }
  
  // Now calculate the sum of the differences squared (Chi2) over all spectra in
  // theSpectra. Display the fit parameters on the standard output and finally,
  // return Chi2 to Minuit.
  double Chi2 = 0.0;
  double Chi2W = 0.0;
  double TheoryPoint;
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    vector<Coordinate> theDiffEXAFS =theSpectra->at(i).SubStructure[0].Spectrum.getDataPoints();
    for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
      TheoryPoint = 0.0;
      for (unsigned int k = 0; k < theSpectra -> at(i).SubStructure.size (); k ++) {
        TheoryPoint += theSpectra -> at (i).SubStructure[k].Spectrum.dChi(theDiffEXAFS[j].x);
      }
      if (theSpectra -> at(i).SubStructure[0].Spectrum.sigmaSqr(j) != 0.0) {
        Chi2 += pow(TheoryPoint - theDiffEXAFS[j].y, 2);
        Chi2W += pow(TheoryPoint - theDiffEXAFS[j].y, 2) 
          // / theSpectra -> at(i).SubStructure[0].Spectrum.sigmaSqr(j);
		/ 1.0e-8;
        //cout << theDiffEXAFS[j].x << "  " << theDiffEXAFS[j].y << "  " << theSpectra -> at(i).SubStructure[0].Spectrum.sigmaSqr(j) << "  " << Chi2W << endl;
      }
    }
  }

  if (VerboseLevel >= MAX_VERBOSE) {
    cout.precision (10);
    cout << Iteration << scientific;  
    for (unsigned int i = 0; i < par.size (); i ++) {
      cout << "  " << par [i];
    }
    cout << "  " << Chi2 << "  " << Chi2W << endl;
  }
  
  return Chi2W;
}

int MgstFcn::setVerbose (int NewLevel) {
  if (NewLevel >= NO_VERBOSE && NewLevel <= MAX_VERBOSE) {
    VerboseLevel = NewLevel;
    return MF_NO_ERROR;
  } else {
    cout << "ERROR: Requested verbose level (" << NewLevel << ") is invalid in MgstFcn::setVerbose" << endl;
    return MF_ERROR;
  }
}
