// MgstSpectrum.h
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
#ifndef MGST_SPECTRUMH
#define MGST_SPECTRUMH

#include "MgstTensor.h"
#include "PathFinder.h"
#include "DiffXasSpectrum.h"

#define MS_NO_ERROR 0
#define MS_ERROR    1

#define DEF_AVE_STEPS 64

// Important: Beamline frame of reference is
//    +z beam propagation direction,
//    +y from hutch floor vertically up towards ceiling
//    +x away from the storage ring in direction orthogonal to +z and +y

typedef struct MS_DeltaS {
  double AveMag1;
  double AveMag2;
  double Ave;
  Path *Chip;
  MS_DeltaS () { AveMag1 = 0.0; AveMag2 = 0.0; Ave = 0.0; Chip = NULL; }
} DeltaS;

typedef struct MS_PrefOrient {
  double x;
  double y;
  double z;
  double Degree;
} PrefOrient;

class MgstSpectrum :
  public MgstTensor, 
  public DiffXasSpectrum,
  public PathFinder {

public:
  MgstSpectrum () { SpectrumReady = false; NumAveSteps = DEF_AVE_STEPS; 
    dChiCheb = gsl_cheb_alloc (CHEB_ORDER); }
  ~MgstSpectrum () { /*gsl_cheb_free (dChiCheb); */}

  // Local declarations of DiffXasSpectrum abstract base functions
  double dChi (double x);
  FFData getDChiPath (int PathNo);
  
  int copySpectrum (MgstSpectrum Input);
  virtual int numSteps (int NewNumSteps) 
    { NumAveSteps = NewNumSteps; SpectrumReady = false; return MS_NO_ERROR; }
  virtual int numSteps () { return NumAveSteps; }
  virtual int addPrefOrientation (double x, double y, double z, double Degree);

  virtual int polarisation (double x, double y, double z) 
    { Polarisation.x = x; Polarisation.y = y;
      Polarisation.z = z; SpectrumReady = false; return MS_NO_ERROR; }
  virtual Coordinate polarisation () { return Polarisation; }
  virtual int magnetisation1 (double x, double y, double z) 
    { Magnetisation1.x = x; Magnetisation1.y = y;
      Magnetisation1.z = z; SpectrumReady = false; return MS_NO_ERROR; }
  virtual Coordinate magnetisation1 () { return Magnetisation1; }
  virtual int magnetisation2 (double x, double y, double z) 
    { Magnetisation2.x = x; Magnetisation2.y = y;
      Magnetisation2.z = z; SpectrumReady = false; return MS_NO_ERROR; }
  virtual Coordinate magnetisation2 () { return Magnetisation2; }

  virtual int mgstCoefficients (double NewLA0, double NewLG1, double NewLE2) { MgstTensor::mgstCoefficients(NewLA0, NewLG1, NewLE2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambdaA0 (double NewLambdaA0) { MgstTensor::lambdaA0 (NewLambdaA0); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambdaG2 (double NewLambdaG2) { MgstTensor::lambdaG2 (NewLambdaG2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambdaE2 (double NewLambdaE2) { MgstTensor::lambdaE2 (NewLambdaE2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambdaD2 (double NewLambdaD2) { MgstTensor::lambdaD2 (NewLambdaD2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda1A0 (double NewLambda1A0) { MgstTensor::lambda1A0 (NewLambda1A0); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda2A0 (double NewLambda2A0) { MgstTensor::lambda2A0 (NewLambda2A0); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda1A2 (double NewLambda1A2) { MgstTensor::lambda1A2 (NewLambda1A2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda2A2 (double NewLambda2A2) { MgstTensor::lambda2A2 (NewLambda2A2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda1G2 (double NewLambda1G2) { MgstTensor::lambda1G2 (NewLambda1G2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda2G2 (double NewLambda2G2) { MgstTensor::lambda2G2 (NewLambda2G2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda3G2 (double NewLambda3G2) { MgstTensor::lambda3G2 (NewLambda3G2); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int lambda4G2 (double NewLambda4G2) { MgstTensor::lambda4G2 (NewLambda4G2); SpectrumReady = false; return MS_NO_ERROR; }    
  virtual int tensorElement (int Element, double NewValue) { MgstTensor::tensorElement (Element, NewValue); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int voigtElement (int Element, double NewValue) { MgstTensor::voigtElement (Element, NewValue); SpectrumReady = false; return MS_NO_ERROR; }
  virtual int halfVoigtElement (int Element, double NewValue) { MgstTensor::halfVoigtElement (Element, NewValue); SpectrumReady = false; return MS_NO_ERROR; }
  int DS (int Path, double NewValue) { if (Path >= 0 && Path < (int)dS.size ()) { dS [Path].Ave = NewValue; SpectrumReady = false; return MS_NO_ERROR; } else { return MS_ERROR; } }

  virtual double lambdaA0 () { return MgstTensor::lambdaA0 (); }
  virtual double lambdaG2 () { return MgstTensor::lambdaG2 (); }
  virtual double lambdaE2 () { return MgstTensor::lambdaE2 (); }
  virtual double lambdaD2 () { return MgstTensor::lambdaD2 (); }
  virtual double lambda1A0 () { return MgstTensor::lambda1A0 (); }
  virtual double lambda2A0 () { return MgstTensor::lambda2A0 (); }
  virtual double lambda1A2 () { return MgstTensor::lambda1A2 (); }
  virtual double lambda2A2 () { return MgstTensor::lambda2A2 (); }
  virtual double lambda1G2 () { return MgstTensor::lambda1G2 (); }
  virtual double lambda2G2 () { return MgstTensor::lambda2G2 (); }
  virtual double lambda3G2 () { return MgstTensor::lambda3G2 (); }
  virtual double lambda4G2 () { return MgstTensor::lambda4G2 (); }
  virtual double tensorElement (int i) { return MgstTensor::tensorElement (i); }
  virtual double voigtElement (int i) { return MgstTensor::voigtElement (i); }
  virtual double halfVoigtElement (int i) {return MgstTensor::halfVoigtElement (i);}
  double DS (int i) { return dS[i].Ave; }
  
  // Get and set functions for the preferential orientations. These are of the
  // same basic template as the path property get and set functions and so can
  // be used by ScriptReader, allowing the preferential orientation properties
  // to be handled by fitparameter and so used as fit parameters.
  virtual int prefX (int Orientation, double NewX);
  virtual int prefY (int Orientation, double NewY);
  virtual int prefZ (int Orientation, double NewZ);
  virtual int prefDeg (int Orientation, double NewDeg);
  virtual double prefX (int Orientation);
  virtual double prefY (int Orientation);
  virtual double prefZ (int Orientation);
  virtual double prefDeg (int Orientation);
  
private:
  void calculateMagnetostriction ();
  void calculateDChi ();
  vector <Coordinate> doDChiCalculation (vector <DeltaS> &dSIn, vector <Coordinate> XCoords);
  Coordinate doDChiCalculation1 (vector <DeltaS> &dSIn, double k);
  static double dChiChebFunction (double k, void *p);
  int rotateVector (Coordinate &theVector, double x, double y, double z);
  void contractPaths ();
  void contractMagnetisation ();
  void contractLegs (vector <Leg*> *PathLegs, double dSMag1, double dSMag2, 
           int theShell, double x, double y, double z, double Pol, 
           Coordinate OldLeg, vector <DeltaS*> *alldS);
  double contractToTensor (Coordinate theLegIn, Coordinate Mag);
  unsigned int NumAveSteps;
  void generateDChiCheb (vector <Coordinate> &Data);
  
  vector <DeltaS> dS;
  vector <PrefOrient> PrefOrientations;
  gsl_cheb_series *dChiCheb;  // A Chebyshev approximation of dChi(k)

  Coordinate Polarisation;    // LINEAR x-ray polarisation wrt beamline frame
  Coordinate Magnetisation1;  // 1st Magnetisation vector wrt beamline frame
  Coordinate Magnetisation2;  // 2nd Magnetisation vector wrt beamline frame
};

#endif // MGST_SPECTRUM_H
