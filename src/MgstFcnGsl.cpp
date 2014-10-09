// MgstFcnGsl.cpp
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
#include <cassert>
#include "Deriv.cpp"

using namespace::std;

// This class and its associated NON-CLASS functions calculate the all the info
// required by GSL's Levenberg-Marquardt fitter, namely the differences vector 
// and the Jacobian matrix. GSL is not quite as sophisticated as Minuit so more
// work has to be done here. Additionally, since GSL is based on C and not C++,
// the functions that are called by GSL may NOT be class members. Instead the 
// void* argument to each function has been used to pass through the class data 
// (in the form of a struct) making them somewhat pseudo-class functions. Much 
// use is made of function pointers here to allow the various spectrum parameter
// functions to be called in an efficient manner.
//

//------------------------------------------------------------------------------
// fitFunction (const gsl_vector *, void *, gsl_vector *) : This function
// calculates the difference vector, which contains the difference between the
// theory and experimental spectra at each of the experiment's data points. The
// spectra are passed into this function via the void *data argument, which is
// of type struct spectra (defined in MgstFcnGsl.h). The difference vector is 
// put in *f, which returns to lower level GSL routines.
//
int fitFunction (const gsl_vector *par, void *data, gsl_vector *f) 
{
//  cout << "diff fn" << endl;
  vector<struct spectrum> *theSpectra = ((vector<struct spectrum> *) data);

  // Use the Index property of each fit parameter in theSpectra to apply the new
  // parameter values to the actal MgstSpectrum objects. 
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    for (unsigned int k = 0; k < theSpectra -> at(i).SubStructure.size (); k ++) {
  
      // Apply the parameters relating to specific scattering paths first...
      for (unsigned int j = 0; j < theSpectra -> at (i).SubStructure[k].Path.size (); j ++) {
        (theSpectra -> at (i).SubStructure[k].Spectrum.*theSpectra -> at (i).SubStructure[k].Path[j].ParamSet)
          (theSpectra -> at (i).SubStructure[k].Path[j].PathNum, 
            gsl_vector_get (par, theSpectra -> at (i).SubStructure[k].Path[j].Index));
//        cout << "Set param: i=" << i << " k=" << k << " j=" << j << " Index=" << theSpectra -> at (i).SubStructure[k].Path[j].Index << " -> " << gsl_vector_get (par, theSpectra -> at (i).SubStructure[k].Path[j].Index) << endl;
      }

      // ... and then the parameter relating to the crystal as a whole
      for (unsigned int j = 0; j < theSpectra -> at (i).SubStructure[k].Crystal.size (); j ++) {
        (theSpectra -> at (i).SubStructure[k].Spectrum.*theSpectra -> at (i).SubStructure[k].Crystal[j].ParamSet)
          (gsl_vector_get (par, theSpectra -> at (i).SubStructure[k].Crystal[j].Index));
//        cout << "Set param: i=" << i << " k=" << k << " j=" << j << " Index=" << theSpectra -> at (i).SubStructure[k].Crystal[j].Index << " -> " << gsl_vector_get (par, theSpectra -> at (i).SubStructure[k].Crystal[j].Index) << endl;
      }
    }
  }

  // Calculate the sum of differences between theory and experiment
  //double Chi = 0.0;
  double ChiW = 0.0;
  double ExptPoint;
  unsigned int Point = 0;
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    vector<Coordinate> theDiffEXAFS =theSpectra->at(i).SubStructure[0].Spectrum.getDataPoints();
    for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
      ExptPoint = 0.0;
      for (unsigned int k = 0; k < theSpectra -> at(i).SubStructure.size (); k ++) {
        ExptPoint += theSpectra -> at (i).SubStructure[k].Spectrum.dChi(theDiffEXAFS[j].x);
      }

      if (theSpectra -> at(i).SubStructure[0].Spectrum.sigmaSqr(j) != 0.0) {
        // Note that the GSL difference vector requires Chi rather than Chi^2 !!!
        //Chi = ExptPoint - theDiffEXAFS[j].y;
        ChiW = (ExptPoint - theDiffEXAFS[j].y)
          / pow (theSpectra -> at(i).SubStructure[0].Spectrum.sigmaSqr(j), 0.5);
//        cout << ExptPoint << "  " << Point << ": i=" << i << " j=" << j << " sig2=" 
//          << theSpectra -> at(i).SubStructure[0].Spectrum.sigmaSqr(j) << " chiW=" << ChiW << endl;
      } else {
        cout << "sigSqr is ZERO in diff fn" << endl;
      }
      gsl_vector_set (f, Point, ChiW);
      Point ++;
    }
  }

  return GSL_SUCCESS;
}


//------------------------------------------------------------------------------
// derivativeFunction (const gsl_vector *, void *, gsl_matrix *) : This function
// calculates the Jacobian matrix, that is a matrix containing the derivative of
// the spectrum with respect to each fit parameter (df/dx_i). The spectra are 
// passed into this function via the void *data argument, which is of type 
// struct spectra (defined in MgstFcnGsl.h).The derivatives are calculated using
// GSL's gsl_deriv_central function. Function pointers to the spectrum's 
// parameter set and get functions are passed through to GSL so it may change
// the parameters (see the dChi function below).
//
int derivativeFunction (const gsl_vector *Arguments, void *data, gsl_matrix *J)
{
//  cout << "deriv fn" << endl;
  vector<struct spectrum> *theSpectra = ((vector<struct spectrum> *) data);
  struct dataCrystal FunctionCrystal;
  struct dataPath FunctionPath;
  gsl_function F;

  double Result, Error;
  double Step = 1.0e-6;
  
  // First initalise the Jacobian (J) so all elements are zero
  for (unsigned int j = 0; j < J -> size1; j ++) {
    for (unsigned int k = 0; k < J -> size2; k ++) {
      gsl_matrix_set (J, j, k, 0.0);
    }
  }
  
  vector<MgstSpectrum> derivSpectra;
  
  // Calculate the derivative of each spectrum i with respect to each parameter 
  // k at each spectrum data point j
  unsigned int Point = 0;
//  unsigned int PointStruct = 0;
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    vector<Coordinate> theDiffEXAFS = theSpectra -> at(i).SubStructure[0].Spectrum.getDataPoints();
    for (unsigned int l = 0; l < theSpectra -> at(i).SubStructure.size(); l ++) {
      FunctionCrystal.theSpectrum.push_back(&theSpectra -> at(i).SubStructure[l].Spectrum);
      FunctionPath.theSpectrum.push_back(&theSpectra -> at(i).SubStructure[l].Spectrum);
    }

    for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
      cout << "\b\b\b\b\b\b\bj = " << j << flush;
      FunctionCrystal.x = theDiffEXAFS[j].x;
      FunctionPath.x = theDiffEXAFS[j].x;

      for (unsigned int l = 0; l < theSpectra -> at(i).SubStructure.size(); l ++) {
        // Now calculate the derivatives for the path dependent parameters
        F.function = &dChiPath;
        FunctionCrystal.CurrentSpectrum = l;
        FunctionPath.CurrentSpectrum = l;
        for (unsigned int k = 0; k < theSpectra -> at(i).SubStructure[l].Path.size (); k ++) {
          FunctionPath.ParamSet = theSpectra -> at(i).SubStructure[l].Path[k].ParamSet;
          FunctionPath.ParamGet = theSpectra -> at(i).SubStructure[l].Path[k].ParamGet;
          FunctionPath.PathNumber = theSpectra -> at(i).SubStructure[l].Path[k].PathNum;
          F.params = &FunctionPath;
          gsl_deriv_central (&F, 0.0, Step, &Result, &Error);
          gsl_matrix_set 
            (J, Point, theSpectra -> at(i).SubStructure[l].Path[k].Index, (gsl_matrix_get (J, Point,  theSpectra -> at(i).SubStructure[l].Path[k].Index) + Result)/Step);
//          cout << "Path   : Point=" << Point << " Index=" << theSpectra -> at(i).SubStructure[l].Path[k].Index << " Result=" << gsl_matrix_get (J, Point,  theSpectra -> at(i).SubStructure[l].Path[k].Index) << "   " << i << "  " << j << "  " << k << "  " << l << endl;
        }

        // And then the derivatives for the crystal dependent parameters
        F.function = &dChiCrystal;
        for (unsigned int k = 0; k < theSpectra -> at(i).SubStructure[l].Crystal.size (); k ++) {
          FunctionCrystal.ParamSet = theSpectra -> at(i).SubStructure[l].Crystal[k].ParamSet;
          FunctionCrystal.ParamGet = theSpectra -> at(i).SubStructure[l].Crystal[k].ParamGet;
          F.params = &FunctionCrystal;
          gsl_deriv_central (&F, 0.0, Step, &Result, &Error);
          gsl_matrix_set 
            (J, Point, theSpectra -> at(i).SubStructure[l].Crystal[k].Index, (gsl_matrix_get (J, Point,  theSpectra -> at(i).SubStructure[l].Crystal[k].Index) + Result/Step));
//          cout << "Crystal: Point=" << Point << " Index=" << theSpectra -> at(i).SubStructure[l].Crystal[k].Index << " Result=" << gsl_matrix_get (J, Point,  theSpectra -> at(i).SubStructure[l].Crystal[k].Index) << "   " << i << "  " << j << "  " << k << "  " << l << endl;
        }
      }
      Point ++;
    }
  }
  cout << "\b\b\b\b\b\b\b" << flush;
  return GSL_SUCCESS;
}




/*int derivativeFunction (const gsl_vector *Arguments, void *data, gsl_matrix *J)
{
  vector<struct spectrum> *theSpectra = ((vector<struct spectrum> *) data);
  crystal_data FunctionCrystal;
  path_data FunctionPath;
  parameter_function F;

  vector<double> Result, Error, x;
  double Step = 1.0;
  
  // First initalise the Jacobian (J) so all elements are zero
  for (unsigned int j = 0; j < J -> size1; j ++) {
    for (unsigned int k = 0; k < J -> size2; k ++) {
      gsl_matrix_set (J, j, k, 0.0);
    }
  }
  
  vector<MgstSpectrum> derivSpectra;

  x.push_back(0.0);
  
  // Calculate the derivative of each spectrum i with respect to each parameter 
  // k at each spectrum data point j
  unsigned int Point = 0;
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    FunctionCrystal.theSpectrum = &theSpectra -> at(i).Spectrum;
    FunctionPath.theSpectrum = &theSpectra -> at(i).Spectrum;
    vector<Coordinate> theDiffEXAFS =theSpectra->at(i).Spectrum.getDataPoints();

    for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
      cout << "\b\b\b\b\b\b\bj = " << j << flush;
      x[0] = theDiffEXAFS[j].x;
      
      // Now calculate the derivatives for the path dependent parameters
      F.function = &dChiPath;
      for (unsigned int k = 0; k < theSpectra -> at(i).Path.size (); k ++) {
        FunctionPath.ParamSet = theSpectra -> at(i).Path[k].ParamSet;
        FunctionPath.ParamGet = theSpectra -> at(i).Path[k].ParamGet;
        FunctionPath.PathNumber = theSpectra -> at(i).Path[k].PathNum;
        F.params = &FunctionPath;
        deriv_parameter_central (&F, x, Step, &Result, &Error);
        gsl_matrix_set 
          (J, Point, theSpectra -> at(i).Path[k].Index, Result[0] / Step);
      }

      // And then the derivatives for the crystal dependent parameters
      F.function = &dChiCrystal;
      for (unsigned int k = 0; k < theSpectra -> at(i).Crystal.size (); k ++) {
        FunctionCrystal.ParamSet = theSpectra -> at(i).Crystal[k].ParamSet;
        FunctionCrystal.ParamGet = theSpectra -> at(i).Crystal[k].ParamGet;
        F.params = &FunctionCrystal;
        deriv_parameter_central (&F, x, Step, &Result, &Error);
        gsl_matrix_set 
          (J, Point, theSpectra ->at(i).Crystal[k].Index, Result[0] / Step);
      }
      Point ++;
    }
  }
  cout << "\b\b\b\b\b\b\b" << flush;
  return GSL_SUCCESS;
}*/


/*int derivativeFunction (const gsl_vector *Arguments, void *data, gsl_matrix *J)
{
  vector<struct spectrum> *theSpectra = ((vector<struct spectrum> *) data);
  crystal_data FunctionCrystal;
  path_data FunctionPath;
  parameter_function F;

  vector<double> Result, Error, x;
  double Step = 1.0;
  
  // First initalise the Jacobian (J) so all elements are zero
  for (unsigned int j = 0; j < J -> size1; j ++) {
    for (unsigned int k = 0; k < J -> size2; k ++) {
      gsl_matrix_set (J, j, k, 0.0);
    }
  }
  
  vector<MgstSpectrum> derivSpectra;
  
  // Calculate the derivative of each spectrum i with respect to each parameter 
  // k at each spectrum data point j
  unsigned int Point = 0;
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    FunctionCrystal.theSpectrum = &theSpectra -> at(i).Spectrum;
    FunctionPath.theSpectrum = &theSpectra -> at(i).Spectrum;
    vector<Coordinate> theDiffEXAFS =theSpectra->at(i).Spectrum.getDataPoints();

    x.clear ();
    for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
      x.push_back(theDiffEXAFS[j].x);
      x.push_back(theDiffEXAFS[j].x);
    }
      
    // Now calculate the derivatives for the path dependent parameters
    F.function = &dChiPath;
    for (unsigned int k = 0; k < theSpectra -> at(i).Path.size (); k ++) {
      FunctionPath.ParamSet = theSpectra -> at(i).Path[k].ParamSet;
      FunctionPath.ParamGet = theSpectra -> at(i).Path[k].ParamGet;
      FunctionPath.PathNumber = theSpectra -> at(i).Path[k].PathNum;
      F.params = &FunctionPath;
      deriv_parameter_central (&F, x, Step, &Result, &Error);

      for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
        gsl_matrix_set 
          (J, Point + j, theSpectra -> at(i).Path[k].Index, Result[j] / Step);
      }
    }

    // And then the derivatives for the crystal dependent parameters
    F.function = &dChiCrystal;
    for (unsigned int k = 0; k < theSpectra -> at(i).Crystal.size (); k ++) {
      FunctionCrystal.ParamSet = theSpectra -> at(i).Crystal[k].ParamSet;
      FunctionCrystal.ParamGet = theSpectra -> at(i).Crystal[k].ParamGet;
      F.params = &FunctionCrystal;
      deriv_parameter_central (&F, x, Step, &Result, &Error);

      for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
        gsl_matrix_set 
          (J, Point + j, theSpectra ->at(i).Crystal[k].Index, Result[j] / Step);
      }
    }
    Point += theDiffEXAFS.size ();
  }
  return GSL_SUCCESS;
}*/


//------------------------------------------------------------------------------
// fitAndDerivativeFunctions (const gsl_vector *, void *, gsl_vector *, 
// gsl_matrix *) : Just calls fitFunction and derivativeFunction.
int fitAndDerivativeFunctions (
  const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J) 
{
  cout << "diff & deriv fn" << endl;
  fitFunction (x, data, f);
  derivativeFunction (x, data, J);
  return GSL_SUCCESS;
}


//------------------------------------------------------------------------------
// dChi (double, void *) : This is the function called by GSL to calculate the
// derivative of the input spectrum with respect to the input parameter. See
// derivateFunction above for more. Basically, the spectrum, and the parameter
// get and set functions are passed in as function pointers so that this 
// function is entirely general whatever the spectrum or parameter function. The
// only requirement is that the parameter get and set functions conform to the
// type found in the struct data definition in MgstFcnGsl.h.
//
double dChiCrystal (double Step, void * data) 
{
  vector<MgstSpectrum*> theSpectra = ((struct dataCrystal *) data) -> theSpectrum;
  unsigned int i = ((struct dataCrystal *) data) -> CurrentSpectrum;
  double x = ((struct dataCrystal *) data) -> x;
  double RtnVal = 0.0;

  (theSpectra[i] ->* ((struct dataCrystal *) data) -> ParamSet) 
    ((theSpectra[i] ->* ((struct dataCrystal *) data) -> ParamGet) () + Step);
  for (unsigned int j = 0; j < theSpectra.size (); j ++) {
    RtnVal += theSpectra[j] -> dChi (x);
  }
  (theSpectra[i] ->* ((struct dataCrystal *) data) -> ParamSet) 
    ((theSpectra[i] ->* ((struct dataCrystal *) data) -> ParamGet) () - Step);
  return RtnVal;
}

double dChiPath (double Step, void * data) 
{
  vector<MgstSpectrum*> theSpectra = ((struct dataPath *) data) -> theSpectrum;
  unsigned int i = ((struct dataPath *) data) -> CurrentSpectrum;
  double x = ((struct dataPath *) data) -> x;
  int PathNumber = ((struct dataPath *) data) -> PathNumber;
  double RtnVal = 0.0;

  (theSpectra[i] ->* ((struct dataPath *) data) -> ParamSet) 
    (PathNumber, (theSpectra[i] ->* 
    ((struct dataPath *) data) -> ParamGet) (PathNumber) + Step);
  for (unsigned int j = 0; j < theSpectra.size (); j ++) {
    RtnVal += theSpectra[j] -> dChi (x);
  }
  (theSpectra[i] ->* ((struct dataPath *) data) -> ParamSet) 
    (PathNumber, (theSpectra[i] ->* 
    ((struct dataPath *) data) -> ParamGet) (PathNumber) - Step);    
  return RtnVal;
}


/*vector<double> dChiCrystal (double Step, vector<double> x, void * data) 
{
  MgstSpectrum *theSpectrum = ((crystal_data *) data) -> theSpectrum;
  vector<double> RtnVals;

  (theSpectrum ->* ((crystal_data *) data) -> ParamSet) 
    ((theSpectrum ->* ((crystal_data *) data) -> ParamGet) () + Step);

  for (unsigned int i = 0; i < x.size (); i ++) {
    RtnVals.push_back (theSpectrum -> dChi (x[i]));
  }
  return RtnVals;
}





vector<double> dChiPath (double Step, vector<double> x, void * data) 
{
  MgstSpectrum *theSpectrum = ((path_data *) data) -> theSpectrum;
  vector<double> RtnVals;
  int PathNumber = ((path_data *) data) -> PathNumber;

  (theSpectrum ->* ((path_data *) data) -> ParamSet) 
    (PathNumber, (theSpectrum ->* 
      ((path_data *) data) -> ParamGet) (PathNumber) + Step);
  
  for (unsigned int i = 0; i < x.size (); i ++) {
    RtnVals.push_back (theSpectrum -> dChi (x[i]));
  }
  return RtnVals;
}*/

//------------------------------------------------------------------------------
// operator() (vector<double> &) : Takes in a vector of doubles containing the
// first guesses to each of the fit parameters, and then uses the GSL non-linear
// scaled Levenberg-Marquardt solver to optimise them.
//
// Only an error code is returned, and the original Guess vector is untouched. 
// The fit results are available in the MgstSpectrum object that was passed by 
// ref to MgstFcn in its constructor. i.e. that spectrum has been optimised.
//
int MgstFcn::operator() (vector<double> &Guess) 
{
  // Prepare the GSL Solver and associated objects. A non-linear solver is used,   // the precise type of which is determined by SOLVER_TYPE, defined in 
  // MgstFcn.h. 
  const size_t NumParameters = Guess.size();

  size_t NumDataPoints = 0;
  double NumIndependentPoints = 0;
  for (unsigned int i = 0; i < theSpectra -> size(); i ++) {
    for (unsigned int k = 0; k < theSpectra -> at(i).SubStructure.size(); k ++) {
      NumDataPoints += theSpectra -> at(i).SubStructure[k].Spectrum.getDataPoints().size();
    }
    NumIndependentPoints += theSpectra -> at(i).SubStructure[0].Spectrum.numIndependentPoints();
  }
  double GuessArr [NumParameters];
  for (unsigned int i = 0; i < NumParameters; i ++) {
    GuessArr[i] = Guess[i];
  }

  const gsl_multifit_fdfsolver_type *SolverType;
  gsl_multifit_fdfsolver *Solver;  
  gsl_multifit_function_fdf FitFunction;
  gsl_matrix *Covariance = gsl_matrix_alloc (NumParameters, NumParameters);
  gsl_vector_view VectorView = gsl_vector_view_array (GuessArr, NumParameters);

  FitFunction.f = &fitFunction;
  FitFunction.df = &derivativeFunction;
  FitFunction.fdf = &fitAndDerivativeFunctions;
  FitFunction.n = NumDataPoints;
  FitFunction.p = NumParameters;
  FitFunction.params = theSpectra;
 
  SolverType = SOLVER_TYPE;
  Solver = gsl_multifit_fdfsolver_alloc(SolverType,NumDataPoints,NumParameters);
  gsl_multifit_fdfsolver_set (Solver, &FitFunction, &VectorView.vector);
  
  // Perform the fitting, one iteration at a time until one of the following
  // conditions is met: The absolute and relative changes in the fit parameters
  // become smaller than SOLVER_TOL, or the max number of allowed iterations,
  // SOLVER_MAX_ITERATIONS, is reached.
  unsigned int Iteration = 0;
  int Status;
  printState (Iteration, NumParameters, Solver);
  do {
    Iteration ++;
    Status = gsl_multifit_fdfsolver_iterate (Solver);
    printState (Iteration, NumParameters, Solver);
    cout << "Status: " << gsl_strerror (Status) << endl;
    if (Status) break;
    Status = gsl_multifit_test_delta (Solver->dx, Solver->x, SOLVER_TOL, SOLVER_TOL);
  } while (Status == GSL_CONTINUE && Iteration < SOLVER_MAX_ITERATIONS 
    && Iteration < MaxIterations);
  cout << "Status: " << gsl_strerror (Status) << endl;


  // This piece of code multiplies the covariance matrix by the fit residuals
  // and divides the result by sigma^2 as explained in section 37.8 of the GSL
  // manual. It is commented out as it is only needed for UN-WEIGHTED least-
  // squares fits. fitFunction() above uses weighted least-squares.
/*  unsigned int Point = 0;
  double NewValue;
  for (unsigned int i = 0; i < theSpectra -> size (); i ++) {
    for (unsigned int l = 0; l < theSpectra -> at(i).SubStructure.size (); l ++) {
      vector<Coordinate> theDiffEXAFS =theSpectra->at(i).SubStructure[l].Spectrum.getDataPoints();
      for (unsigned int j = 0; j < theDiffEXAFS.size (); j ++) {
        for (unsigned int k = 0; k < NumParameters; k ++) {
          NewValue = gsl_matrix_get (Solver -> J, Point, k);
//          cout << "Point=" << Point << "  k=" << k << "  l=" << l << "  " << NewValue << flush;
          NewValue *= gsl_vector_get (Solver -> f, Point);
          NewValue /= pow (theSpectra -> at(i).SubStructure[l].Spectrum.sigmaSqr(j), 0.5);
          gsl_matrix_set (Solver -> J, Point, k, NewValue);
//          cout << "  " << NewValue << endl;
        }
        Point ++;
      }
    }
  }*/

  

  // Output all the fit parameters with their associated error.
  gsl_multifit_covar (Solver -> J, 0.0, Covariance);
#define FIT(i) gsl_vector_get (Solver -> x, i)
#define ERR(i) sqrt (gsl_matrix_get (Covariance, i, i))
  cout << "Covariance Matrix:" << endl;
  for (unsigned int i = 0; i < NumParameters; i ++) {
    for (unsigned int j = 0; j < NumParameters; j ++) {
      cout << gsl_matrix_get (Covariance, i, j) << ", " << flush;
    }
    cout << endl;
  }
  
  {
    double chi = gsl_blas_dnrm2 (Solver -> f);
    double dof = NumIndependentPoints - double(NumParameters);
    double c = GSL_MAX_DBL (1, chi / sqrt (dof));
    
    cout << "Chi^2/DOF = " << pow(chi, 2) / dof << endl;
    cout << "DOF = " << dof << "  c = " << c << endl;
    for (unsigned int j = 0; j < NumParameters; j ++) {
      cout << "Parameter " << j << " = " << FIT(j) << " +/- "<< c*ERR(j)<< endl;
    }
  }
  cout << "Status: " << gsl_strerror (Status) << endl;
 
  // Clean up the memory and exit
  gsl_multifit_fdfsolver_free (Solver);
  gsl_matrix_free (Covariance);
  return 0;
}


//------------------------------------------------------------------------------
// printState (int, int, gsl_multifit_fdfsolver) : Prints to standard output the
// fit parameters and chi value for the current optimisation iteration.
//
void MgstFcn::printState (int i, int p, gsl_multifit_fdfsolver *Solver) 
{
  cout << "Iteration: " << i << "  x = ";
  for (int j = 0; j < p; j ++) {
    cout << gsl_vector_get (Solver -> x, j) << ", ";
  }
  cout << "\b\b   Chi = " << gsl_blas_dnrm2 (Solver -> f) << endl;
}














