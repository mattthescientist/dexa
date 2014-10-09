// DiffXasSpectrum.cpp
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

// This is an abstract base class that provides the low-level I/O functionality
// required to construct a DiffXAS spectrum. The exact type of DiffXAS spectrum
// (magnetostrictive, electostrictive, thermal,etc.) is determined by inheriting
// this class and defining the 'dChi' function, which should implement the
// appropriate form of the DiffXAS equation for the system under study. Other
// definitions must include the two forms of 'chi', and 'getDChiPath'.
//
// This file is divided into four sections. The first for routines relating to
// Experimental spectra, the second to Theoretical data files from FEFF, the
// third for GSL Chebyshev routines used for modelling the theoretical phase, 
// Phi(k), and amplitude, A(k), functions, and the fourth for Fourier filtering
// routines. See the headers for each of these sections for more information on
// the routines they contain.
//
// It is important to note that the () operator is overloaded in this class to
// provide a short-hand link to dChi().
//

#include "DiffXasSpectrum.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include "Resample.cpp"
#include "Fourier.h"

// GSL includes for Chebyshev polynomials
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

// GSL includes for Fast-Fourier transforms (mixed-radix)
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>

// Define PI if it hasn't already been given
#ifndef PI
  #define PI (4.0 * atan (1.0))
#endif

// Define the filename sufficies used when saving F filter intermediate data
#define PATH_K_SUFFIX "_KSpace"
#define PATH_R_SUFFIX "_RSpace"
#define PATH_Q_SUFFIX "_QSpace"

#define HANN_NAME_K2R "HannWindowKtoR.dat"
#define HANN_NAME_R2Q "HannWindowRtoQ.dat"

using namespace::std;

//------------------------------------------------------------------------------
// Default constructor : Just set default values for class variables.
//
DiffXasSpectrum::DiffXasSpectrum () {
  thedK = 0.0;
  thedE = 0.0;
  Verbose = DX_MIN_VERBOSE;

  // Disable explicit definition of q-space fit range
  fqMin = DX_DISABLED; 
  fqMax = DX_DISABLED;
  
  // Disable Fourier filter
  UseFourierFilter = false;
}


//------------------------------------------------------------------------------
// Default destructor : Free all GSL Chebyshevs before destroying the class
//
DiffXasSpectrum::~DiffXasSpectrum () {
  for (unsigned int i = 0; i < ScatteringPaths.size (); i ++) {
    gsl_cheb_free (ScatteringPaths[i].ACheb);
    gsl_cheb_free (ScatteringPaths[i].PhiCheb);
    gsl_cheb_free (ScatteringPaths[i].PCheb);
    ScatteringPaths.clear ();
  }
}


//******************************************************************************
// EXPERIMENT SECTION : Load spectrum, Save spectrum
//******************************************************************************
// This section simply provides routines to load experimental DiffXAS
// spectra taken on ID24.'loadSpectrum' handles files in a free format. The user
// must specify which columns contain the X and Y data. Columns are delimited
// according to standard C++ white-space rules. The user must also declare how
// big the file header is i.e. how many lines must be discarded before the first
// data line is reached. The data must then follow in one chunk with no blank 
// lines or additional comments in the middle. Equally, there must be no
// comments after the data. If any of these rules are not met, just edit your
// data file before running it through here.
//

//------------------------------------------------------------------------------
// loadSpectrum (char*, int, int, int) : Loads a spectrum from the data file 
// specified by arg1. arg2 determines how many lines to discard at the beginning
// of the file, and arg3 and arg4 the columns for X and Y data respectively.
//
int DiffXasSpectrum::loadSpectrum (
  const char *Filename, int Header, int XCol, int YCol) {

  int i;
  string theString;
  Coordinate NextPoint;
  vector<double> NextLine;
  vector<Coordinate> Residuals;
  double NextCol;
  istringstream iss;
  
  vector<double> XCoords, YCoords, ResidualsX, ResidualsY;

  // Open the spectrum data file. Return an error if it doesn't open.
  ifstream SpectrumFile (Filename, ios::in);
  if (! SpectrumFile.is_open()) {
    return DX_ERROR;
  }
  
  // Ditch the header
  for (i = 0; i < Header; i ++) {
    getline (SpectrumFile, theString);
  }

  // Read the data and store the required components for X and Y
  // The data is assummed to be in k-space.
  DataPoints.clear ();
  while (!SpectrumFile.eof ()) {
    getline (SpectrumFile, theString);
    iss.str (theString);
    if (theString[0] != '\0') {         // Catch empty lines eg. just before EOF
      while (!iss.eof ()) {
        iss >> NextCol;
        NextLine.push_back (NextCol);
      }
      iss.clear ();
      XCoords.push_back (NextLine[XCol - 1]);
      YCoords.push_back (NextLine[YCol - 1]);
      ResidualsX.push_back (NextLine [XCol - 1]);
      NextLine.clear ();
    }
  }
  SpectrumFile.close ();
  
  // Now, given that DiffXAS on ID24 is taken in constant steps in E-space, 
  // resample the spectrum so that it has constant steps in k-space. This is
  // required for the FFT. See the description of resample() for details on how
  // this is done. The spacing between data points is the final argument of
  // resample () and is selected so that the maximum signal frequency in the
  // resampled spectrum is that corresponding to an effective photoelctron
  // scattering radius of MAX_RADIUS.
  if (Verbose >= DX_MIN_VERBOSE)
    cout << "...Resampling spectrum to constant k steps" << endl;
  int Err = resample (&XCoords, &YCoords, &ResidualsY, 2.0/(PI*MAX_RADIUS));
  if (Err != DX_NO_ERROR) { return Err; }
  for (unsigned int i = 0; i < XCoords.size (); i ++) {
    NextPoint.x = XCoords[i];
    NextPoint.y = YCoords[i];
    RawData.push_back (NextPoint);
    AllDataPoints.push_back (NextPoint);
    DataPoints.push_back (NextPoint);
  }

  // The residuals returned by resample() contain the statistical noise on the
  // spectrum. They now need to be processed so that they may be used to weight
  // the chi^2 in the fitting algorithm. This processing requires three steps:
  //
  //   1. Take the abs value of the residuals to get the basic noise amplitude.
  //
  //   2. Assuming the noise is white, perform a bandwidth correction to account
  //      for the noise at radii less than MAX_RADIUS, which wasn't returned by
  //      by resample() in the residuals.
  //
  //        Noise *=       R_nyquist
  //                 ----------------------
  //                 R_nyquist - MAX_RADIUS
  //
  //      where the R_nyquist is the effective photoelectron scattering radius
  //      obtained from the spectrum's Nyquist frequency. Since the density of
  //      data points increases with k, the highest sampling rate in k-space will
  //      be that between the last two data points.
  //
  //   3. Correct for intrinsic noise averaging that takes place when 
  //      resampling the spectrum into constant steps in k-space. In essance,
  //      this means dividing by the number of raw data points that were
  //      averaged in order to obtain each resampled point. This is a fn of k:
  //
  //        n(k) =  ds    =  /2 * ALPHA * k\ ds  =    4 * ALPHA * k
  //               -----    |--------------|       --------------------
  //               dk(k)    \     dE      /        PI * dE * MAX_RADIUS
  //
  //      where ds is the resampling step size and dk is the spacing between data
  //      points at k. dK may be expressed in terms of dE, which is a constant.
  if (Verbose >= DX_MIN_VERBOSE)
    cout << "...Determining statistical noise in data" << endl;
  double dKMin = (XCoords[XCoords.size() - 1] - XCoords[XCoords.size() - 2]);
  double RNyquist = (PI * 0.5) / dKMin;
  double dE = ALPHA * (pow(XCoords[1], 2) - pow(XCoords[0], 2));
  for (unsigned int i = 0; i < ResidualsX.size (); i ++) {
    NextPoint.x = ResidualsX[i];
    NextPoint.y = ResidualsY[i];
    Residuals.push_back (NextPoint);
    /* 1 */  ResidualsY[i] = abs(ResidualsY[i]);
    /* 2 */  ResidualsY[i] *= RNyquist / (RNyquist - MAX_RADIUS);
//    /* 3 */  ResidualsY[i] /= (4.0 * ALPHA * ResidualsX[i]) / (PI*dE*MAX_RADIUS);
  }
  if (Verbose >= DX_MAX_VERBOSE) {
    cout << "...Saving resampling debug information" << endl;
    cout << "    RNyquist = " << RNyquist 
         << ", MAX_RADIUS = " << MAX_RADIUS << ", dE = " << dE << endl;
    saveData ("ExptResampledToK.dat", AllDataPoints);
    saveData ("ExptResampledResiduals.dat", Residuals);
  }
  
  //  Finally, take the processed residuals and convert them to a sigma^2 array
  //  that can be used to weight chi^2 in the fitting algorithm.
  if (Verbose >= DX_MIN_VERBOSE)
    cout << "...Calculating sigma^2 for weighted least-squares" << endl;
  processNoise (&ResidualsX, &ResidualsY, MAX_RADIUS, MAX_RADIUS, &XCoords);
  SigmaSqr.clear ();
  for (unsigned int i = 0; i < ResidualsX.size (); i ++) {
    NextPoint.x = ResidualsX[i];
    NextPoint.y = pow (ResidualsY[i], 2.0);
    SigmaSqr.push_back (NextPoint);
    AllSigmaSqr.push_back (NextPoint);
  }
  if (Verbose >= DX_MAX_VERBOSE) {
    cout << "...Saving sigma sqr debug information" << endl;
    saveData ("ExptResampledSigmaSqr.dat", SigmaSqr);
  }
  return DX_NO_ERROR;
}


//------------------------------------------------------------------------------
// saveSpectrum (char*, int, int, int) : Saves the experimental spectrum. This
// can be useful if some processing has been performed on the data since loading
// it i.e. Fourier filtering. This function is really just a wrapper for
// saveData(), needed because DataPoints are not accessible outside the class.
//
int DiffXasSpectrum::saveSpectrum (const char *Filename) {
  return saveData (Filename, DataPoints);
}


//******************************************************************************
// THEORY SECTION : Path manipulation, dChi calculation, file i/o
//******************************************************************************
// This section provides the functions necessary to load scattering path data in
// from FEFF's chipNNNN.dat files. Paths may be added with the either of the two 
// forms of addPath, modified with updatePath, or removed with deletePath. These
// provide all the necessary functionality to construct a spectrum from theory.
// However, in order to obtain specific information on individual scattering
// paths, FEFF's feff.inp and paths.dat must be processed. Operations relating
// to this are found in PathFinder.cpp.
//

//------------------------------------------------------------------------------
// addPath (char*, int) : Loads the feffnnnn.dat file specified by arg1. arg2 
// determines how many lines to discard at the beginning of the file. This 
// should be the lines up to but NOT INCLUDING the line containing "nleg, deg, 
// reff, rnrmav(bohr), edge". This usually means that 8 lines must be discarded.
//
int DiffXasSpectrum::addPath (const char *Filename, int Header) {

  ifstream DataFile (Filename, ios::in);
  if (! DataFile.is_open()) {
    return DX_ERROR;
  }
  
  // Ditch the first part of the header
  string theString;
  for (int i = 0; i < Header; i ++) {
    getline (DataFile, theString);
  }
  
  // Make sure that we've been given a 'feff' theory data file
  theString.assign (Filename, 0, 4); 
  if (theString == "feff") { readFeffFile (DataFile, Filename); }
  else {
    cout << "Error: Unrecognised file type in addPath.\n"
      << "Filename must start with \"feff\"." << endl;
    DataFile.close ();
    return DX_ERROR;
  }
  DataFile.close ();
  return DX_NO_ERROR;
}


//------------------------------------------------------------------------------
// readFeffFile (ifstream &, const char *, string) : Given an input file stream,
// loaded with one of FEFF's feffnnnn.dat files, this funtion will extract all
// the information required to construct the absorption fine-struction function,
// chi. Most parameters may be directly read from the file, but the amplitude
// and phase functions, A(k) and Phi(k) respectively, must be calculated from
// their constituent factors:
//
//   A_j(k) = N_j F_j(k) exp[-2R_j / lambda(k)]     Phi_j(k) = p(2phc) + p(feff)
//            ---------------------------------
//                          kR_j^2
//
// For more details, see "EXAFS analysis using FEFF and FEFFIT" by Matthew
// Newville, J. Synchrotron Rad. 8 (2001) 96-100.
//
// At the end of the function, Chebyshev polynomials of A(k), Phi(k), and P(k)
// are generated. These are used for all subsequent evaluations of each of these
// functions are are necessary for two reasons. Firstly, they allow evaluation
// at an arbirary k, and more importantly, easily accommodate changes in the 
// edge position as is necessary when applying dE or dk corrections.
//
int DiffXasSpectrum::readFeffFile(ifstream &FeffFile, const char *FeffFilename){
  string theString;
  Path NewPath;
  vector<double> NextLine;
  double NextCol;
  istringstream iss;
  double Degeneracy, Radius;
  
  // Get the path index by finding the the 'nnnn' part of 'feffnnnn.dat'
  theString.assign (FeffFilename, 4, 4); 
  NewPath.Index = atoi (theString.c_str ());
  
  // Read the necessary path properties from the feffnnnn.dat header. These are
  // N_j, R_j, S02_j, and sigma_j^2. The name of the feffnnnn.dat file is also
  // stored so that any scattering path found by PathFinder may be associated
  // with its appropriate feffnnnn.dat information.
  string Temp;
  getline (FeffFile, theString); iss.str (theString);
  iss >> Temp >> Degeneracy >> Radius;
  do {
    getline (FeffFile, theString); iss.str (theString);
    iss >> Temp;
  } while (Temp != "k" && !FeffFile.eof ());
  if (FeffFile.eof ()) { 
    cout << "Error reading " << FeffFilename << "." << endl;
    return DX_ERROR;
  }
  NewPath.S = Radius * 2.0;
  NewPath.S02 = 1.0;
  NewPath.ChipName = FeffFilename;
  NewPath.Sig2 = 0.0;
  
  // Read the A(k) and Phi(k) functions and store them in NewPath. Each must be
  // calculated from other factors in the feffnnnn.dat file. See the function
  // header above for more details.
  getline (FeffFile, theString);
  while (!FeffFile.eof ()) {
    getline (FeffFile, theString);
    iss.str (theString);
    if (theString[0] != '\0') {   // Check the line isn't blank e.g. near EOF
      NextLine.clear ();
      while (!iss.eof ()) {
        iss >> NextCol;
        NextLine.push_back (NextCol);
      }
      iss.clear ();
      
      NewPath.A.push_back (
        Degeneracy * NextLine[FEFF_MAG_COL-1] * NextLine[FEFF_RED_FAC_COL-1] 
        * exp (-2.0 *Radius/NextLine[FEFF_LAMBDA_COL-1])/(NextLine[FEFF_K_COL-1]
        * Radius * Radius));
      NewPath.Phi.push_back (
        NextLine[FEFF_REAL_2PHC_COL-1] + NextLine[FEFF_PHASE_COL-1]);
      NewPath.K.push_back (NextLine[FEFF_K_COL-1]);
      NewPath.P.push_back (NextLine[FEFF_REAL_P_COL -1]);
    }
  }

  // Finally, generate Chebyshev approximations of each of the A(k), Phi(k), and   // P(k) functions, and store the complete path in the ScatteringPaths class
  // vector.
  NewPath.ACheb = gsl_cheb_alloc (CHEB_ORDER);
  NewPath.PhiCheb = gsl_cheb_alloc (CHEB_ORDER);
  NewPath.PCheb = gsl_cheb_alloc (CHEB_ORDER);
  evalTheoryChebyshevs (&NewPath);
  
  ScatteringPaths.push_back (NewPath);
  return DX_NO_ERROR;
}


//------------------------------------------------------------------------------
// deletePath (char*) : Deletes the path with chipnnnn.dat filename matching the
// input parameter. GSL Chebyshevs are freed before deleting the path itself.

int DiffXasSpectrum::deletePath (char *ChipName) {
  for (unsigned int i = 0; i < ScatteringPaths.size (); i ++) {
    if (ScatteringPaths[i].ChipName == ChipName) {
      gsl_cheb_free (ScatteringPaths[i].ACheb);
      gsl_cheb_free (ScatteringPaths[i].PhiCheb);
      gsl_cheb_free (ScatteringPaths[i].PCheb);
      ScatteringPaths.erase (ScatteringPaths.begin () + i);
      return DX_NO_ERROR;
    }
  }
  return DX_ERROR;
} 


//------------------------------------------------------------------------------
// updatePath (string, string, double) : Finds the stored path matching ChipName
// and replaces the property specified by Field with NewValue. Field must be
// S02, Sig2, S, R, or dR, and are case-insensitive.
//
int DiffXasSpectrum::updatePath (string ChipName, string Field, double NewValue) {
  for (unsigned int i = 0; i < ScatteringPaths.size (); i ++) {
    if (ScatteringPaths[i].ChipName == ChipName) {
      if (Field == "s02" || Field == "S02") {
        ScatteringPaths[i].S02 = NewValue;
        return DX_NO_ERROR;
      } else if (Field == "sig2" || Field == "Sig2" || Field == "SIG2" 
              || Field == "SiG2" || Field == "SIg2" || Field == "sIG2") {
        ScatteringPaths[i].Sig2 = NewValue;
        return DX_NO_ERROR;
      } else if (Field == "s" || Field == "S") {
        ScatteringPaths[i].S = NewValue;
        return DX_NO_ERROR;
      } else if (Field == "r" || Field == "R") {
        ScatteringPaths[i].S = NewValue * 2.0;
        return DX_NO_ERROR;
      } else if (Field == "dR" || Field == "dr" || Field == "Dr" ||Field=="DR"){
        ScatteringPaths[i].dR = NewValue;
        return DX_NO_ERROR;
      } else {
        return DX_ERROR;
      }
    }
  }
  return DX_ERROR;
}


//------------------------------------------------------------------------------
// getPath (int) : Finds and returns the path with index matching arg1. This is
// a protected function, which is not called within the DiffXasSpectrum class, 
// but which is available to be called by any child class. It is NOT public 
// since the path is returned by reference rather than by value, making it 
// possible for an external function to damage DiffXasSpectrum's internal path 
// list. External functions should therefore call getPaths() or getPathj(), 
// which return by value.
//
Path* DiffXasSpectrum::getPath (int Index) {
  for (unsigned int i = 0; i < ScatteringPaths.size (); i ++) {
    if (ScatteringPaths[i].Index == Index) return &ScatteringPaths[i];
  }
  return NULL;
}


//------------------------------------------------------------------------------
// operator () : Overload the () operator so that it calls dChi.
//
double DiffXasSpectrum::operator() (double kIn) {
  return dChi (kIn);
}

//------------------------------------------------------------------------------
// saveTheory (char*) : Saves the calculated dChi spectrum. Data points are
// calculated and saved at the same positions in k space as those loaded from
// the experimental spectrum. Column 1 is for the k ordinate, 2 for Chi, and
// the following columns for contributions from each scattering path. Data is
// saved in scientific notation to 5 decimal places, and columns are delimited
// with a single tab character.
//
int DiffXasSpectrum::saveTheory (const char *Filename, vector<Coordinate> &DiffXas) {

  vector < vector <Coordinate> > AllKData;
  vector < vector <Coordinate> > AllRData;
  vector < vector <Coordinate> > AllQData;
  FFData PathData = getDChiPath (0);
  
  AllKData.push_back (getTheoryKData ());
  AllRData.push_back (getTheoryRData ());
  AllQData.push_back (getTheoryQData ());
  DiffXas = AllQData[0];
  for (unsigned int j = 0; j < ScatteringPaths.size (); j ++) {
    PathData = getDChiPath (j);
    AllKData.push_back (PathData.KData);
    AllRData.push_back (PathData.RData);
    AllQData.push_back (PathData.QData);
  }

  ostringstream oss;  
  if (UseFourierFilter && Verbose >= DX_MIN_VERBOSE) {
    // Save the K-space data
    if (Verbose >= DX_MED_VERBOSE) {
      oss << Filename << PATH_K_SUFFIX;
      writeTheoryData (oss.str(), AllKData);
    }

    // Save the R-space data
    oss.clear (); oss.str ("");
    oss << Filename << PATH_R_SUFFIX;
    writeTheoryData (oss.str(), AllRData);
    
    // Save the Q-space data
//    oss.clear (); oss.str ("");
//    oss << Filename << PATH_Q_SUFFIX;
//    writeTheoryData (oss.str(), AllQData);
  } else {
    // If we're not using the Fourier Filter, just save the QData without suffix
    writeTheoryData (Filename, AllQData);
  }
  return DX_NO_ERROR;
}


int DiffXasSpectrum::writeTheoryData (string Filename, vector < vector <Coordinate> > &Data) {

  unsigned int i, j;
  ofstream SpectrumFile (Filename.c_str(), ios::out);
  if (! SpectrumFile.is_open()) { return DX_ERROR; }
  for (i = 0; i < Data[0].size (); i ++) {
    SpectrumFile << showpos << scientific << Data[0][i].x;
    for (j = 0; j < Data.size (); j ++) {
      SpectrumFile << showpos << scientific << '\t' << Data[j][i].y;
    }
    SpectrumFile << endl;
  }
  SpectrumFile.close ();
  return DX_NO_ERROR;
}


//------------------------------------------------------------------------------
// chi (double) : Calculates the standard EXAFS fine-structure function. chi 
// (double, int) is called to obtain the individual contributions from each
// scattering path, which are then summed here to give the total chi.
//
double DiffXasSpectrum::chi (double kIn) {

  double Chi = 0.0;
  for (unsigned int j = 0; j < ScatteringPaths.size (); j ++) {
    Chi += chi (kIn, j);
  }
  return Chi;
}


//------------------------------------------------------------------------------
// chi (double, int) : Calculates the standard EXAFS fine-structure contribution
// at k Ang^-1 for scattering path PathNo. The calculation is based upon eqn 3
// in "EXAFS analysis using FEFF and FEFFIT" by Matthew Newville, J. Synchrotron
// Rad. 8 (2001) 96-100, but where contributions from the third and fourth
// cumulants have been neglected.
//
double DiffXasSpectrum::chi (double kIn, int PathNo) {
  double chi = 0.0;
  double A, P,Phi, Sig2, S, k, dR;

  k = pow(((kIn * kIn * 3.81) - dE()) / 3.81, 0.5);
  k = k + dK();
  
  int j = PathNo;
  A    = gsl_cheb_eval (ScatteringPaths[j].ACheb, k) * ScatteringPaths[j].S02;
  P    = gsl_cheb_eval (ScatteringPaths[j].PCheb, k);
  Phi  = gsl_cheb_eval (ScatteringPaths[j].PhiCheb, k);
  Sig2 = ScatteringPaths[j].Sig2;
  S    = ScatteringPaths[j].S;
  dR   = ScatteringPaths[j].dR;

  A = A * pow (S * 0.5, 2) / pow ((S * 0.5) + dR, 2);

  chi += A * sin((S * k) + Phi + (2.0 * P * dR) - (8.0 * P * Sig2 / S))
    * exp (-2.0 * P * P * Sig2);
  return chi;
}

  
//******************************************************************************
// CHEBYSHEV SECTION
//******************************************************************************
// The functions here use the GSL library to generate Chebyshev polynomial
// approximations of the scattering phase and amplitude functions loaded from
// FEFF's chipNNNN.dat files. This is done so that, the theoretical chi function
// may be accurately evaluated at the experimental data points, which are almost
// certainly differ from those generated by FEFF.
//

//------------------------------------------------------------------------------
// evalChebyshevs (Path*) : Evaluates the Chebyshev approximations for A_j(k)
// and Phi_j(k) for the Path specified at arg1. These Chebyshevs are stored as
// part of the path data itself.
//
void DiffXasSpectrum::evalTheoryChebyshevs (Path *thePath) {
  gsl_function F;

  F.params = thePath;  
  F.function = &DiffXasSpectrum::amplitudeFunction;
  gsl_cheb_init (thePath->ACheb, &F, thePath->K[0], 
    thePath->K[thePath->K.size()-1]);

  F.function = &DiffXasSpectrum::phaseFunction;
  gsl_cheb_init (thePath->PhiCheb, &F, thePath->K[0], 
    thePath->K[thePath->K.size()-1]);
    
  if (thePath->P.size() > 0) {
    F.function = &DiffXasSpectrum::pFunction;
    gsl_cheb_init (thePath->PCheb, &F, thePath->K[0], 
      thePath->K[thePath->K.size()-1]);
  }
}


//------------------------------------------------------------------------------
// amplitudeFunction (double, void*) : Used by GSL's gsl_cheb_init to evaluate
// A_j(k) at any value k between the 0.0 and CHIP_MAX_K. At this stage, points 
// between those specified in the chipNNNN.dat files are determined by linear 
// interpolation.
//
double DiffXasSpectrum::amplitudeFunction (double k, void *p) {
  int UpperPoint, LowerPoint, MidPoint;
  
  Path *thePath = (Path *)p;
  LowerPoint = 0; UpperPoint = thePath -> K.size () - 1;
  
  if (k >= thePath->K[UpperPoint]) { return thePath -> A[UpperPoint]; } 
  else if (k <= thePath->K[LowerPoint]) { return thePath -> A[LowerPoint]; }

  do {
    MidPoint = int((UpperPoint - LowerPoint)/2) + LowerPoint;
    if (MidPoint == UpperPoint || MidPoint == LowerPoint) {
      return thePath->A[LowerPoint] + 
        ((k - thePath -> K[LowerPoint]) /
        ((thePath -> K[UpperPoint] - thePath -> K[LowerPoint]))) *
        (thePath->A[UpperPoint] - thePath->A[LowerPoint]);
    }

    if (k > thePath -> K[MidPoint]) {
      LowerPoint = MidPoint;
    } else {
      UpperPoint = MidPoint;
    }
  } while (true);
}


//------------------------------------------------------------------------------
// phaseFunction (double, void*) : Used by GSL's gsl_cheb_init to evaluate
// Phi_j(k) at any value k between the 0.0 and CHIP_MAX_K. At this stage,
// points between those specified in the chipNNNN.dat files are determined by
// linear interpolation.
//
double DiffXasSpectrum::phaseFunction (double k, void *p) {
  int UpperPoint, LowerPoint, MidPoint;
  
  Path *thePath = (Path *)p;
  LowerPoint = 0; UpperPoint = thePath -> K.size () - 1;

  if (k >= thePath -> K[UpperPoint]) { return thePath -> Phi[UpperPoint]; } 
  else if (k <= thePath -> K[LowerPoint]) { return thePath -> Phi[LowerPoint]; }

  do {
    MidPoint = int((UpperPoint - LowerPoint)/2) + LowerPoint;
    if (MidPoint == UpperPoint || MidPoint == LowerPoint) {
      return thePath->Phi[LowerPoint] + 
        ((k - thePath -> K[LowerPoint]) /
        ((thePath -> K[UpperPoint] - thePath -> K[LowerPoint]))) *
        (thePath->Phi[UpperPoint] - thePath->Phi[LowerPoint]);
    }

    if (k > thePath -> K[MidPoint]) {
      LowerPoint = MidPoint;
    } else {
      UpperPoint = MidPoint;
    }
  } while (true);
}


//------------------------------------------------------------------------------
// pFunction (double, void*) : Used by GSL's gsl_cheb_init to evaluate
// P_j(k) at any value k between the 0.0 and CHIP_MAX_K. At this stage,
// points between those specified in the chipNNNN.dat files are determined by
// linear interpolation.
//
double DiffXasSpectrum::pFunction (double k, void *p) {
  int UpperPoint, LowerPoint, MidPoint;
  
  Path *thePath = (Path *)p;
  LowerPoint = 0; UpperPoint = thePath -> K.size () - 1;

  if (k >= thePath -> K[UpperPoint]) { return thePath -> P[UpperPoint]; } 
  else if (k <= thePath -> K[LowerPoint]) { return thePath -> P[LowerPoint]; }

  do {
    MidPoint = int((UpperPoint - LowerPoint) / 2) + LowerPoint;
    if (MidPoint == UpperPoint || MidPoint == LowerPoint) {
      return thePath->P[LowerPoint] + 
        ((k - thePath -> K[LowerPoint]) /
        ((thePath -> K[UpperPoint] - thePath -> K[LowerPoint]))) *
        (thePath->P[UpperPoint] - thePath->P[LowerPoint]);
    }

    if (k > thePath -> K[MidPoint]) {
      LowerPoint = MidPoint;
    } else {
      UpperPoint = MidPoint;
    }
  } while (true);
}


//******************************************************************************
// FOURIER FILTER SECTION
//******************************************************************************


double DiffXasSpectrum::numIndependentPoints () {
  double KRange = 0.0;
  double RRange = 0.0;
  for (unsigned int i = 0; i < fkMin.size (); i ++) {
    KRange += fkMax[i] - fkMin[i];
    RRange += frMax[i] - frMin[i];
  }
  if (!UseFourierFilter || KRange == 0.0 || RRange == 0.0) {
    return DataPoints.size ();
  } else {
    return 2.0 * KRange * RRange / PI;
  }
}


//------------------------------------------------------------------------------
// setFilter (double, double, double, double, double, double) : Sets the Fourier
// transform's Hann window properties. The first three arguments specify the
// properties of the forward transform from k- to R-space, and the last three,
// the back transform from R- to q-space. For more details on how they are used
// see the comments associated with hannWindow() below.
//
void DiffXasSpectrum::setFilter (double nkMin, double nkMax, double ndK, 
  double nrMin, double nrMax, double ndR) {
  fkMin.push_back (nkMin); fkMax.push_back (nkMax); fdK.push_back (ndK);
  frMin.push_back (nrMin); frMax.push_back (nrMax); fdR.push_back (ndR);
  UseFourierFilter = true;
  
/*  vector<double> XNoise, YNoise, XSample;
  for (unsigned int i = 0; i < Residuals.size (); i ++) {
    XNoise.push_back (Residuals[i].x);
    YNoise.push_back (Residuals[i].y);
  }
  for (unsigned int i = 0; i < DataPoints.size (); i ++) {
    XSample.push_back (RawData[i].x);
  }
  processNoise (&XNoise, &YNoise, frMax, MAX_RADIUS, &XSample);
  vector <Coordinate> AllSigmaSqr;  
  SigmaSqr.clear ();
  Coordinate NextPoint;
  for (unsigned int i = 0; i < XNoise.size (); i ++) {
    NextPoint.x = XNoise[i];
    NextPoint.y = pow (YNoise[i], 2.0);
    AllSigmaSqr.push_back (NextPoint);
  }*/
  
  AllDataPoints = RawData;
  fourierFilter (fkMin, fkMax, fdK, frMin, frMax, fdR, AllDataPoints, ExptFFData);
/*  if (fqMin != DX_DISABLED && fqMax != DX_DISABLED) {
    setFitRange (fqMin, fqMax);
  } else {
    setFitRange (fkMin, fkMax);    
  }*/
  DataPoints.clear ();
  SigmaSqr.clear ();
  for (unsigned int i = 0; i < AllDataPoints.size (); i ++) {
    for (unsigned int j = 0; j < fkMin.size (); j ++) {
      if (AllDataPoints[i].x >= fkMin[j] && AllDataPoints[i].x <= fkMax[j]) {
        DataPoints.push_back (AllDataPoints[i]);
        SigmaSqr.push_back (AllSigmaSqr[i]);
      }
    }
  } 
}


void DiffXasSpectrum::setFitRange (double nqMin, double nqMax) {
  fqMin = nqMin;
  fqMax = nqMax;

  DataPoints.clear ();
  for (unsigned int i = 0; i < AllDataPoints.size (); i ++) {
    if (AllDataPoints[i].x >= fqMin && AllDataPoints[i].x <= fqMax) {
      DataPoints.push_back (AllDataPoints[i]);
    }
  }
}

//------------------------------------------------------------------------------
// fourierFilter (double, double, double, double, double, double, vector 
// <Coordinate> &, FFData &) : Extracts a specific range of scattering path 
// lengths from the DiffXAS spectrum for further analysis. The spectrum is first
// transformed from k- to R-space using data in the range kMin <= k <= kMax. 
// This data is obtained by applying a Hann window of width parameter dK. 
// Likewise, the R-space data is then backtransformed to q-space in the range 
// rMin <=R <=rMax. Again, a Hann window is used; this time with width parameter
// dR. The original input spectrum is OVERWRITTEN with the filtered data.
//
// All the intermediate spectra (windowed and un-windowed k-, R-, and q-space)
// are returned in the FFData object at the final argument.
//
int DiffXasSpectrum::fourierFilter (vector<double> kMin, vector<double> kMax, 
  vector<double> dK, vector<double> rMin, vector<double> rMax, 
  vector<double> dR, vector <Coordinate> &Data, FFData &Output) {

  // If the Filter isn't being used, just fill Output with the contents of Data
  // and return.
  if (!UseFourierFilter) { 
    for (unsigned int i = 0; i < Data.size (); i ++) {
      Output.KData.push_back (Data[i]);
      Output.QData.push_back (Data[i]);
      Output.KWindowed.push_back (Data[i]);
    }
    return DX_NO_ERROR; 
  }
  
  Fourier FFT;
  vector <Coordinate> TotalFT;

  // First the forward FT from k- to R-space
  Output.KData = Data;
  if (Verbose >= DX_MED_VERBOSE) {
    FFT.hannWindow (kMin, kMax, dK, Data, HANN_NAME_K2R);
  } else {
    FFT.hannWindow (kMin, kMax, dK, Data);
  }
  Output.KWindowed = Data;
  FFT.transform (Data, FORWARD, Data[0].x, Output.RData);
  
  // Then the backward FT from R- to q-space
  Output.RWindowed = Output.RData;
  if (Verbose >= DX_MED_VERBOSE) {
    FFT.hannWindowComplex (rMin, rMax, dR, Output.RWindowed, HANN_NAME_R2Q);
  } else {
    FFT.hannWindowComplex (rMin, rMax, dR, Output.RWindowed);
  }
  FFT.transform (Output.RWindowed, BACKWARD, Data[0].x, Output.QData);
  
  // Overwrite the input k-space spectrum with the filtered q-space spectrum
  Data.clear ();
  for (unsigned int i = 0; i < Output.QData.size(); i ++) {
    if (Output.QData[i].x >= EXAFS_MIN_K && Output.QData[i].x <= EXAFS_MAX_K) {
      Data.push_back (Output.QData[i]);
    }
  }
  Output.QData = Data;
  return DX_NO_ERROR;
}


//------------------------------------------------------------------------------
// Fourier filter GET functions. Pretty self-explanatory except for the R-space
// spectra. For each case (Expt and Theory spectra), get...RData constructs the
// total FT from the complex class data, get...RDataReal returns just the real
// component of the FT, and get...RDataComplex just the complex component.
vector <Coordinate> DiffXasSpectrum::getExptKData () {
  return ExptFFData.KData;
}

vector <Coordinate> DiffXasSpectrum::getExptRData () {
  vector <Coordinate> TotalFt;
  Fourier FFT;
  FFT.getTotalFT (ExptFFData.RData, TotalFt);
  return TotalFt;
}

vector <Coordinate> DiffXasSpectrum::getExptRDataReal () {
  vector <Coordinate> RealFTComponent;
  for (unsigned int i = 0; i < ExptFFData.RData.size() / 2; i ++) {
    RealFTComponent.push_back (ExptFFData.RData [i]);
  }
  return RealFTComponent;
}

vector <Coordinate> DiffXasSpectrum::getExptRDataComplex () {
  vector <Coordinate> ComplexFTComponent;
  for (unsigned int i = ExptFFData.RData.size() / 2; i < ExptFFData.RData.size(); i ++) {
    ComplexFTComponent.push_back (ExptFFData.RData [i]);
  }
  return ComplexFTComponent;
}

vector <Coordinate> DiffXasSpectrum::getExptQData () {
  return ExptFFData.QData;
}

vector <Coordinate> DiffXasSpectrum::getExptRWindowed () {
  vector <Coordinate> TotalFt;
  Fourier FFT;
  FFT.getTotalFT (ExptFFData.RWindowed, TotalFt);
  return TotalFt;
}

vector <Coordinate> DiffXasSpectrum::getTheoryKData () {
  return TheoryFFData.KData;
}

double DiffXasSpectrum::getTheoryKData (int DataPoint) {
  return TheoryFFData.KData[DataPoint].y;
}

vector <Coordinate> DiffXasSpectrum::getTheoryRData () {
  vector <Coordinate> TotalFt;
  Fourier FFT;
  FFT.getTotalFT (TheoryFFData.RData, TotalFt);
  return TotalFt;
}

vector <Coordinate> DiffXasSpectrum::getTheoryRDataReal () {
  vector <Coordinate> RealFTComponent;
  for (unsigned int i = 0; i < TheoryFFData.RData.size() / 2; i ++) {
    RealFTComponent.push_back (TheoryFFData.RData [i]);
  }
  return RealFTComponent;
}

vector <Coordinate> DiffXasSpectrum::getTheoryRDataComplex () {
  vector <Coordinate> ComplexFTComponent;
  for (unsigned int i = TheoryFFData.RData.size() / 2; i < TheoryFFData.RData.size(); i ++) {
    ComplexFTComponent.push_back (TheoryFFData.RData [i]);
  }
  return ComplexFTComponent;
}


vector <Coordinate> DiffXasSpectrum::getTheoryQData () {
  return TheoryFFData.QData;
}

vector <Coordinate> DiffXasSpectrum::getTheoryRWindowed () {
/*  vector <Coordinate> TotalFt;
  Fourier FFT;
  FFT.getTotalFT (TheoryFFData.RWindowed, TotalFt); 
  return TotalFt;*/
  return TheoryFFData.RWindowed;
}


//------------------------------------------------------------------------------
// saveSpectrum (char*, int, int, int) : saves the EXPERIMENTAL data stored
// within the DiffXasSpectrum object. This is useful for when the experimental
// data has been Fourier filtered and thus differs from the data in the spectrum
// loaded at the beginning of an analysis.
//
int DiffXasSpectrum::saveData (const char *Filename, vector<Coordinate> Data) {

  ofstream SpectrumFile (Filename, ios::out);
  if (! SpectrumFile.is_open()) {
    return DX_ERROR;
  }
  
  SpectrumFile.precision (5);
  for (unsigned int i = 0; i < Data.size (); i ++) {
    SpectrumFile << showpos << scientific 
      << Data[i].x << '\t' << Data[i].y << '\t' << endl;
  }
  SpectrumFile.close();
  return DX_NO_ERROR;
}


//------------------------------------------------------------------------------
// loadSigmaSqr (const char *, int, int, int) :
//
int DiffXasSpectrum::loadSigmaSqr (
  const char *Filename, int Header, int XCol, int YCol) {
  
  int i;
  ifstream SpectrumFile (Filename, ios::in);
  string theString;
  Coordinate NextPoint;
  vector<double> NextLine;
  vector<Coordinate> LoadedSigmaSqr;
  double NextCol;
  istringstream iss;
  
  vector<Coordinate> SigmaSqrData;
  Coordinate NewSigmaSqr;

  if (! SpectrumFile.is_open()) {
    return DX_ERROR;
  }
  
  // Ditch the header
  for (i = 0; i < Header; i ++) {
    getline (SpectrumFile, theString);
  }

  // Read the data and store the required components for X and Y
  // The data is assummed to be in k-space.
  while (!SpectrumFile.eof ()) {
    getline (SpectrumFile, theString);
    iss.str (theString);
    if (theString[0] != '\0') {         // Catch empty lines eg. just before EOF
      while (!iss.eof ()) {
        iss >> NextCol;
        NextLine.push_back (NextCol);
      }
      iss.clear ();
      NewSigmaSqr.x = NextLine[XCol - 1];
      NewSigmaSqr.y = NextLine[YCol - 1];
      SigmaSqrData.push_back (NewSigmaSqr);
      NextLine.clear ();
    }
  }
  SpectrumFile.close ();
  
  // Now create a Chebyshev representation of the sigma sqr to allow accurate
  // interpolation at the DiffXAS data points.
  gsl_cheb_series *SigmaSqrCheb;
  gsl_function F;
  
  SigmaSqrCheb = gsl_cheb_alloc (CHEB_ORDER);
  F.params = &SigmaSqrData;
  F.function = &DiffXasSpectrum::SigmaSqrChebFunction;
  gsl_cheb_init (
    SigmaSqrCheb, &F, SigmaSqrData[0].x, SigmaSqrData[SigmaSqrData.size()-1].x);

  // Finally, evaluate sigma sqr at each of the DiffXAS data points and replace
  // the SigmaSqr class array.
  for (unsigned int j = 0; j < SigmaSqr.size(); j ++) {
    NewSigmaSqr.x = SigmaSqr[j].x;
    NewSigmaSqr.y = abs(gsl_cheb_eval (SigmaSqrCheb, NewSigmaSqr.x));
//    cout << NewSigmaSqr.x << "  " << NewSigmaSqr.y << endl;
    LoadedSigmaSqr.push_back (NewSigmaSqr);
  }
  SigmaSqr.clear();
  AllSigmaSqr.clear();
  SigmaSqr = LoadedSigmaSqr;
  AllSigmaSqr = LoadedSigmaSqr;
    
  return DX_NO_ERROR;
}


//------------------------------------------------------------------------------
// SigmaSqrChebFunction (double, void *) : 
//
double DiffXasSpectrum::SigmaSqrChebFunction (double k, void *p) {
  int UpperPoint, LowerPoint, MidPoint;
  
  vector <Coordinate> *dChiPts = (vector <Coordinate> *)p;
  LowerPoint = 0; UpperPoint = dChiPts -> size () - 1;
  
  if (k >= dChiPts->at(UpperPoint).x) { return dChiPts->at(UpperPoint).y; } 
  else if (k <= dChiPts->at(LowerPoint).x) { return dChiPts->at(LowerPoint).y; }

  do {
    MidPoint = int((UpperPoint - LowerPoint)/2) + LowerPoint;
    if (MidPoint == UpperPoint || MidPoint == LowerPoint) {
      return dChiPts->at(LowerPoint).y + 
        ((k - dChiPts->at(LowerPoint).x) /
        ((dChiPts->at(UpperPoint).x - dChiPts->at(LowerPoint).x))) *
        (dChiPts->at(UpperPoint).y - dChiPts->at(LowerPoint).y);
    }

    if (k > dChiPts->at(MidPoint).x) {
      LowerPoint = MidPoint;
    } else {
      UpperPoint = MidPoint;
    }
  } while (true);
}
