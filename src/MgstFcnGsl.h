// MgstFcnGsl.h
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
#ifndef MF_MGST_FIT_GSL_H
#define MF_MGST_FIT_GSL_H

#define DEF_MAX_ITERATIONS 250

#include "MgstSpectrum.h"
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

#ifndef NULL
  #define NULL 0
#endif

#define SOLVER_TYPE gsl_multifit_fdfsolver_lmsder
#define SOLVER_TOL 1.0e-8
#define SOLVER_MAX_ITERATIONS 500

class MgstFcn {

public:

  MgstFcn (std::vector <struct spectrum> *s) 
    {theSpectra = s; MaxIterations = DEF_MAX_ITERATIONS; }
  ~MgstFcn () {};
  
  int operator() (std::vector<double> &Guess);
  int iterations (int iter) { MaxIterations = iter; return 0; }

private:
  unsigned int MaxIterations;
  void printState (int i, int p, gsl_multifit_fdfsolver *Solver);
  vector <struct spectrum> *theSpectra;
};

struct dataCrystal { 
  vector<MgstSpectrum*> theSpectrum; 
  int (MgstSpectrum::*ParamSet)(double);
  double (MgstSpectrum::*ParamGet)();
  double x;
  unsigned int CurrentSpectrum;
  int MSIndex;
};

struct dataPath { 
  vector<MgstSpectrum*> theSpectrum; 
  int (MgstSpectrum::*ParamSet)(int, double);
  double (MgstSpectrum::*ParamGet)(int);
  double x;
  unsigned int CurrentSpectrum;
  int PathNumber;
};


/*struct crystal_data_struct { 
  MgstSpectrum *theSpectrum; 
  int (MgstSpectrum::*ParamSet)(double);
  double (MgstSpectrum::*ParamGet)();
};

struct path_data_struct { 
  MgstSpectrum *theSpectrum; 
  int (MgstSpectrum::*ParamSet)(int, double);
  double (MgstSpectrum::*ParamGet)(int);
  int PathNumber;
};*/

typedef struct path_data_struct path_data;
typedef struct crystal_data_struct crystal_data;

int fitFunction (const gsl_vector *x, void *data, gsl_vector *f);
int derivativeFunction (const gsl_vector *x, void *data, gsl_matrix *J);
int fitAndDerivativeFunctions 
  (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J);
double dChiCrystal (double Step, void * data);
double dChiPath (double Step, void * data);
vector<double> dChiCrystal (double Step, /*vector<*/double/*>*/ x, void * data);
vector<double> dChiPath (double Step, /*vector<*/double/*>*/ x, void * data);

#endif // MF_MGST_FIT_GSL_H
