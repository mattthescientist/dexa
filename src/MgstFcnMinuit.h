// MgstFcnMinuit.h
// Copyright (C) 2007-2009, 2014 Matthew Ruffoni
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
#ifndef MF_MGST_FCN_H
#define MF_MGST_FCN_H

#define DEF_MAX_ITERATIONS 10000

#define MF_NO_ERROR 0
#define MF_ERROR 1

#include <vector>
#include <iostream>
#include "ThermalSpectrum.h"
#include "Minuit/FCNBase.h"
#include "ScriptReader.h"

using namespace::std;

class MgstFcn : public FCNBase {

public:
  MgstFcn (vector <struct spectrum> *Spectra) {
    theSpectra = Spectra; MaxIterations = DEF_MAX_ITERATIONS; 
    VerboseLevel = MIN_VERBOSE; theErrorDef = 0.0;
  }
  ~MgstFcn () {}
  
  virtual double up() const { return theErrorDef; }
  virtual double operator() (const vector<double>&) const;
  
  int iterations (int iter) { MaxIterations = iter; return 0; }
  void setErrorDef (double NewDef) { theErrorDef = NewDef; }
  int setVerbose (int NewLevel);

private:
  int MaxIterations;
  double theErrorDef;
  int VerboseLevel;
  vector <struct spectrum> *theSpectra;
};

#endif // MF_MGST_FCN_H
