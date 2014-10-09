// DiffXasSpectrum.h
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
#ifndef DX_DIFFXAS_SPECTRUM_H
#define DX_DIFFXAS_SPECTRUM_H

// DiffXasSpectrum error codes
#define DX_NO_ERROR 0
#define DX_ERROR    1
#define DX_DISABLED -1.0

// Verbose Levels
#define DX_NO_VERBOSE  0
#define DX_MIN_VERBOSE 1
#define DX_MED_VERBOSE 2
#define DX_MAX_VERBOSE 3

// Properties of FEFF's feffNNNN.dat files
#define FEFF_K_COL         1
#define FEFF_REAL_2PHC_COL 2
#define FEFF_MAG_COL       3
#define FEFF_PHASE_COL     4
#define FEFF_RED_FAC_COL   5
#define FEFF_LAMBDA_COL    6
#define FEFF_REAL_P_COL    7

#define EXAFS_MIN_K   0.1  // Angstroms^-1
#define EXAFS_MAX_K   20.0 // Angstroms^-1
#define EXAFS_K_STEP  0.05 // Angstroms^-1

// GSL Chebyshev constants
#define CHEB_ORDER    200
#define CHEB_BG_ORDER 3

#define MAX_RADIUS 17 // Angstroms
#define ALPHA      3.81

// Debug Output Filenames
//#define FT_R_SPACE 

#include <vector>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#include <iostream>

using namespace::std;

//------------------------------------------------------------------------------
// TypeDefs
// dSj is not used with all types of DiffXAS signal. For instance, in the
// case of magnetostriction, dSj for each path is calculated from the
// magnetostriction tensor.

typedef struct DX_ChipFile {
  std::string ChipName;     // Name of the chipnnnn.dat file for this path
  std::vector <double> A;   // Aj(k)    Amplitude function
  std::vector <double> Phi; // Phij(k)  Phase function
  std::vector <double> K;   // k        Photoelectron wavevector
  std::vector <double> P;   // p        p = real[p] + i/lambda (in feffnnnn.dat)
  std::vector <double> Lam; // Lambda
  double S02;               // S02      Shake-off
  double Sig2;              // Sigma^2  Debye-Waller factor
  double S;                 // Sj       Path length. Equals 2Rj
  double dR;                // dR       Change in path's R found from EXAFS fit
  double dSig2;
  int Index;

  gsl_cheb_series *ACheb;   // A Chebyshev approximation of A_j(k)
  gsl_cheb_series *PhiCheb; // A Chebyshev approximation of Phi_j(k)
  gsl_cheb_series *PCheb;   // A Chebyshev approximation of P_j(k)
  gsl_cheb_series *LamCheb; // A Chebyshev approximation of Lambda_j(k)

  DX_ChipFile () {
	  S02 = 1.0;
	  Sig2 = 0.0;
	  S = 0.0;
	  dR = 0.0;
	  dSig2 = 0.0;
	  Index = 0;
	  ACheb = NULL;
	  PhiCheb = NULL;
	  PCheb = NULL;
  	  LamCheb = NULL;
  }
} Path;

typedef struct DX_Coordinate {
  double x;
  double y;
  double z;
  DX_Coordinate () { x = 0.0; y = 0.0; z = 0.0; }
} Coordinate;

// DX_FFData the intermediate Fourier filter spectra that are generated by
// fourierFilter(). This allows the calling function to persistantly store this
// information experimental and theory spectra.
typedef struct DX_FFData {
  vector <Coordinate> KData;          // K-space spectrum. Real data storage
  vector <Coordinate> RData;          // R-space spectrum. Complex data storage
  vector <Coordinate> QData;          // Q-space spectrum. Real data storage
  vector <Coordinate> KWindowed;      // K-space after Hann windowing. Real
  vector <Coordinate> RWindowed;      // R-space after Hann windowing. Complex
} FFData;

//------------------------------------------------------------------------------
// MgstSpectrum Class

class DiffXasSpectrum {

public:
  DiffXasSpectrum ();          // Default constructor just defines BG Chebyshev
  virtual ~DiffXasSpectrum (); // Default destructor must free GSL Chebyshevs
  
  // Public Commands
  int loadSpectrum (const char *Filename, int Header, int XCol, int YCol);
  int addPath (const char *ChipFilename, int Header);
  int deletePath (char *Chipfile);
  int updatePath (string ChipName, string Field, double NewValue);
  int saveTheory (const char *Filename, vector<Coordinate> &DiffXas);
  int saveSpectrum (const char *Filename);
  double operator() (double x);

  // Abstract functions that must be implemented in a child class
  virtual double dChi (double x) = 0;
  virtual double chi (double x);
  virtual double chi (double x, int PathNo);
  virtual FFData getDChiPath (int PathNo) = 0;

  // Public GET Routines
  vector<Coordinate> getDataPoints ()  { return DataPoints; }
  vector<Coordinate> getSigmaSqr ()    { return SigmaSqr; }
  Path               getPathj (int j)  { return ScatteringPaths[j]; }
  int                numPaths ()       { return ScatteringPaths.size (); }
  virtual double     dK ()             { return thedK; }
  virtual double     dE ()             { return thedE; }
  virtual double     dSig2 (int Shell) { return ScatteringPaths[Shell].dSig2; }
  double numIndependentPoints ();
  double sigmaSqr (int Point) { return SigmaSqr[Point].y; }
  
  // GET routines for the Fourier Filter information
  vector <Coordinate> getExptKData ();
  vector <Coordinate> getExptRData ();
  vector <Coordinate> getExptQData ();
  vector <Coordinate> getExptRDataReal ();
  vector <Coordinate> getExptRDataComplex ();
  vector <Coordinate> getExptKWindowed () { return ExptFFData.KWindowed; }
  vector <Coordinate> getExptRWindowed ();
  virtual vector <Coordinate> getTheoryKData ();
  virtual double getTheoryKData (int DataPoint);
  vector <Coordinate> getTheoryRData ();
  vector <Coordinate> getTheoryQData ();
  vector <Coordinate> getTheoryRDataReal ();
  vector <Coordinate> getTheoryRDataComplex ();
  vector <Coordinate> getTheoryKWindowed () { return TheoryFFData.KWindowed; }
  vector <Coordinate> getTheoryRWindowed ();
  bool usingFourierFilter () { return UseFourierFilter; }

  // Public SET Routines  
  virtual void dK (double a) { thedK = a; SpectrumReady = false; }
  virtual int dE (double a)  { thedE = a; SpectrumReady = false; return 0; }
  virtual int dSig2 (int Shell, double NewdSig2);
  void setFilter 
    (double kMin, double kMax, double dK, double rMin, double rMax, double dR);
  void setFitRange (double nqMin, double nqMax);
  int loadSigmaSqr (const char* Filename, int Header, int XCol, int YCol);
  void setVerbose (int NewLevel) { Verbose = NewLevel; }
  
 
protected:
  Path* getPath (int Index);
  vector <Path> ScatteringPaths;      // Theory from FEFF's chipnnnn.dat files
  vector<double> fkMin, fkMax, fdK, frMin, frMax, fdR;  // Fourier filter parameters
  double fqMin, fqMax;  // Fit range for back-transformed (q-space) data
  int fourierFilter (vector <double> kMin, vector <double> kMax, vector <double> dK,
    vector <double> rMin, vector <double> rMax, vector <double> dR, vector <Coordinate> &Data,FFData &Data2);
  bool SpectrumReady;  
  bool UseFourierFilter;
  int saveData (const char *Filename, vector<Coordinate> Data);

  // Fourier Filter data
  FFData ExptFFData;                 // Expt Fourier filter diagnostic spectra
  FFData TheoryFFData;               // Theory Fourier filter diagnostic spectra
  vector <Coordinate> RawData;       // Raw experimental DiffEXAFS data
  vector <Coordinate> AllDataPoints; // All Fourier filtered points
  vector <Coordinate> DataPoints;    // F filtered data within given range
  vector <Coordinate> SigmaSqr;      // Residuals processed into sig^2 noise
  vector <Coordinate> AllSigmaSqr;   // Sigma sqr for all data points

private:
  double thedK;                      // Global edge correction parameters
  double thedE;
  int Verbose;
  
  // Chebyshev evaluation functions  
  void evalTheoryChebyshevs (Path *thePath);
  void evalExptChebyshev ();
  static double amplitudeFunction (double k, void *params);
  static double phaseFunction (double k, void *params);
  static double pFunction (double k, void *params);
  static double lambdaFunction (double k, void *params);
  static double exptFunction (double k, void *params);
  static double SigmaSqrChebFunction (double k, void *p);
  int writeTheoryData (string Filename, vector < vector <Coordinate> > &Data);
  
  int readFeffFile (ifstream &FeffFile, const char *FeffFilename);
};

#endif // DX_DIFFXAS_SPECTRUM_H
