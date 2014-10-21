// ThermalSpectrum.cpp
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
#include "ThermalSpectrum.h"

#include <cmath>

// This class inherits and so completes the DiffXasSpectrum abstract base class
// by providing 'dChi' to calculate DiffXAS signals for THERMAL MODULATION.
//

//------------------------------------------------------------------------------
// calculateDChi () : Calculates the entire DiffEXAFS dChi function between
// EXAFS_MIN_K <= k <= EXAFS_MAX_K Ang^-1. Once calculated, the spectrum is
// Fourier filtered and a Chebyshev approximation of the function generated and
// stored for later use by dChi(double) below. Since the Fourier filter must act
// on an entire spectrum, and with data points spaced evenly in k, a stored
// Chebyshev is needed to allow dChi to be obtained for any arbitrary k, such as
// those at which experimental data points exist.
//
void ThermalSpectrum::calculateDChi () {

  vector <Coordinate> temp;
  vector <Coordinate> dChi = doDChiCalculation (ScatteringPaths, temp);

  // The entire dChi(k) spectrum now exists in dChi. Now do the post-processing.
  // First, in order to obtain FT R-space data on the same scale as the
  // experimental data, it is necessary to resample dChi so that the data points
  // are co-incident with the experimental points. This is done with the aid of
  // a Chebyshev polynomial representation of dChi.
  generateDChiCheb (dChi);
  vector <Coordinate> dChiSampled;
  Coordinate NewPoint;
  for (unsigned int i = 0; i < RawData.size(); i ++) {
    NewPoint.x = RawData[i].x;
    NewPoint.y = gsl_cheb_eval (dChiCheb, NewPoint.x);
    dChiSampled.push_back (NewPoint);
  }
  if (dChiSampled.size () == 0) { dChiSampled = dChi; }

  // Now Fourier filter the resampled dChi and regenerate the the dChi Chebyshev
  // so that it contains the filtered spectrum. Since it is dChiSampled that is
  // filtered, the information in TheoryFFData will be comparable with that in
  // ExptFFData.
  fourierFilter (fkMin, fkMax, fdK, frMin, frMax, fdR, dChiSampled, TheoryFFData);
  generateDChiCheb (dChiSampled);
}


//------------------------------------------------------------------------------
// doDChiCalculation (...) : This function is a wrapper for doDchiCalculation1
// below. It allows dChi to be calculated at either regular intervals of
// EXAFS_K_STEP between EXAFS_K_MIN and EXAFS_K_MAX, or at specific points. If
// XCoords is blank at the input, it will be the former, if not, it will be the
// latter, where the points will be taken from the XCoords.x property.
//
vector <Coordinate> ThermalSpectrum::doDChiCalculation (vector <Path> paths, vector <Coordinate> XCoords) {
  vector <Coordinate> dChi;

  double k, kIn;
  if (XCoords.size () == 0) {
    for (kIn = EXAFS_MIN_K; k < EXAFS_MAX_K; kIn += EXAFS_K_STEP) {
      // Only proceed with the calculation if k > 0.0 after any edge shifts
      if (((kIn * kIn * ALPHA) - dE()) / ALPHA - pow(dK(),2) > 0.0) {
        k = pow(((kIn * kIn * ALPHA) - dE()) / ALPHA, 0.5) - dK();

        // Calculate dChi here. Pass in kIn, rather than k, as that must be
        // stored in the x ordinate of the output data points.
        dChi.push_back (doDChiCalculation1 (paths, kIn));
      }
    }
  } else {
    for (unsigned int j = 0; j < XCoords.size (); j ++) {
      k = XCoords[j].x;
      dChi.push_back (doDChiCalculation1 (paths, k));
    }
  }
  return dChi;
}


//------------------------------------------------------------------------------
// doDChiCalculation1 (vector <DeltaS> &) : This is the function that actually
// does the dChi calculations. It is called by doDChiCalculation above.
// The input vector contains all the scattering legs to be considered for the
// calculation in hand. In the case of calculateDChi, this is all the legs. For
// getDChiPath, it will only be the legs pertaining to one particular path.
//
// Note two things: Firstly, the calculation of dChi works with the half path
// length Rj rather than the full length Sj, despite the fact that all the
// magnetostrictive strains are calculated for the latter. This is simply to
// make it easier to implement the EXAFS function from Eq.(3) of M. Newville,
// J. Synchrotron Radiat. 8, 96, which gives Chi(k) for Rj. Rj is simply (1/2)Sj
//
// Secondly, the calculation of dChi does not implement the Taylor series form
// of the DiffXAS structure function here, but instead explicitly calculates the
// difference in Chi as chi'(k) - chi''(k). Since it is easy to obtain both chi'
// and chi'' it is better to work this way as it avoids the limitations of the
// first order Taylor expansion. See lab book 20/10/2008 for more details.
//
Coordinate ThermalSpectrum::doDChiCalculation1 (vector <Path> paths, double kIn){

  double k;
  Coordinate NewDChi;

  // Calculate the DiffXAS at each k by summing contributions over each path j
  // Take the edge shift parameters into account when calculating the value of
  // k to use for the calculation. Take special care to check that any dE does
  // not result in us calculating dChi for k <= 0.0
  NewDChi.y = 0.0;
  if (((kIn * kIn * ALPHA) - dE()) / ALPHA - pow(dK(),2) > 0.0) {
    k = pow(((kIn * kIn * ALPHA) - dE()) / ALPHA, 0.5) - dK();
    for (unsigned int j = 0; j < paths.size (); j ++) {

      //----------------------------------------------------------------------
      // Standard EXAFS components. These are the data read from the feffnnnn
      // files. See DiffXasSpectrum::readFeffFile(...) for more information.
      //
      double A    = gsl_cheb_eval (paths[j].ACheb, k) * paths[j].S02;
      double P    = gsl_cheb_eval (paths[j].PCheb, k);
      double Phi  = gsl_cheb_eval (paths[j].PhiCheb, k);
      double Lam  = gsl_cheb_eval (paths[j].LamCheb, k);
      double Sig2 = paths[j].Sig2;
      double S    = paths[j].S;

      //----------------------------------------------------------------------
      // DiffEXAFS components
      //
      // Rj is half the photoelectron scattering path length S. In the case of
      // single scattering paths, this is the atomic coordination shell radius
      //
      // dRj is the modification to Rj that was found from a conventional XAS
      // analysis. i.e. is the difference between the 'real' Rj and that used
      // by FEFF when calculating the data in the feffnnnn.dat files.
      //
      // dSj is the change in scattering path length due to thermal strain.
      // Note that dS[j].Ave is the NORMALISED change in path length. i.e. it
      // is the strain. It is therefore multiplied by the path length to get
      // the actual atomic displacement.
      //
      // Cj is a COMMON amplitude factor, which contains all the components of
      // the EXAFS function that do not change between measurements the
      // different thermal measurements at T' and T''.
      //
      // Dj is the amplitude factor that DEPENDS upon the sample temperature,
      // and so changes between the two measurements at T' and T''.
      //
      // Qj is the phase component of EXAFS signal. Like Dj, it is dependent
      // upon the thermal state, and so changes between the two measurements.
      //
      // Each variable is assigned the parameters for T', the first
      // thermal state, upon declaration.
      //
      double Rj  = S * 0.5;
      double dRj = paths[j].dR;
      double dSj = paths[j].dS * (S + 2.0 * dRj);
      double dSig2j = paths[j].dSig2;
      double Dj  =
		  A * pow (Rj, 2) * exp (
			  -2.0 * P * P * Sig2
			  + 2.0 * Sig2 / (Lam * Lam)
			  - 2.0 * dRj / Lam
			  + 4.0 * Sig2 / (Rj * Lam)
		  ) * pow (dRj + Rj, -2.0);
      double Qj  =
		  (2.0 * k * Rj)
		  + Phi
		  + (2.0 * P * dRj)
		  - (4 * P * Sig2 / Rj)
		  - (4.0 * P * Sig2 / Lam);

      // Now change the thermal state to T'' and recalculate the dependent
      // components of the DiffEXAFS function (i.e. Dj and Qj)
      dRj += dSj / 2;
      Sig2 += dSig2j;
      double Dj2 =
		  A * pow (Rj, 2) * exp (
			  -2.0 * P * P * Sig2
			  + 2.0 * Sig2 / (Lam * Lam)
			  - 2.0 * dRj / Lam
			  + 4.0 * Sig2 / (Rj * Lam)
		  ) * pow (dRj + Rj, -2.0);
      double Qj2 =
		  (2.0 * k * Rj)
		  + Phi
		  + (2.0 * P * dRj)
		  - (4 * P * Sig2 / Rj)
          - (4.0 * P * Sig2 / Lam);

      // Finally, take the difference between the two measurements to get the
      // overall DiffEXAFS signal.
      NewDChi.y +=  ( (Dj * sin (Qj)) - (Dj2 * sin (Qj2)) );
    }
  }
  NewDChi.x = kIn;
  return NewDChi;
}


//------------------------------------------------------------------------------
// generateDChiCheb (vector <Coordinate>) : Generates a Chebyshev polynomial
// representation of the input data and stores the result in the class variable
// dChiCheb.
//
void ThermalSpectrum::generateDChiCheb (vector <Coordinate> &Data) {
  gsl_cheb_free (dChiCheb);
  dChiCheb = gsl_cheb_alloc (CHEB_ORDER);
  gsl_function F;
  F.params = &Data;
  F.function = &ThermalSpectrum::dChiChebFunction;
  gsl_cheb_init (dChiCheb, &F, Data[0].x, Data[Data.size() - 1].x);
}


//------------------------------------------------------------------------------
// dChiChebFunction (double, void *) : Returns dChi(k) for any arbitrary k by
// using linear interpolation. This function is used by GSL when generating a
// Chebyshev approximation of dChi. Thereafter, the Chebyshev is used to obtain
// dChi(k) rather than this function.
//
double ThermalSpectrum::dChiChebFunction (double k, void *p) {
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
  return 0.0; // Never runs, but stops compiler warnings
}


//------------------------------------------------------------------------------
// dChi (double) : returns dChi(k) by evaluating the Chebyshev approximation of
// dChi at k.
double ThermalSpectrum::dChi (double k) {

  // First make sure that all the strains have been calculated. If
  // they haven't, or need updating, call 'contractPaths' to do that.
  if (!SpectrumReady) {
    calculateDChi ();
    SpectrumReady = true;
  }
  // Then obtain dChi(k) from the dChiCheb Chebyshev
  if (((k * k * 3.81) - dE()) / 3.81 - pow(dK(),2) <= 0.0) return 0.0;
  return gsl_cheb_eval (dChiCheb, k);
}


//------------------------------------------------------------------------------
// getDChiPath (int) : returns dChi_j(k). That is, the contribution to the total
// dChi signal by path j only. This can be used for diagnostic/debug purposes to
// see the effect of individual paths on the overall dChi. As with
// calculateDChi() above, the actual calculation is done in the sub-function
// doDChiCalculation() below. The difference is that only the scattering legs of
// path j are passed to it rather than all the legs.
//
FFData ThermalSpectrum::getDChiPath (int PathNo) {

  // First make sure that all the strains have been calculated. If
  // they haven't, or need updating, call 'contractPaths' to do that.
  if (!SpectrumReady) {
    calculateDChi ();
    SpectrumReady = true;
  }

  // Since we are only returning the dChi component for a single scattering
  // path, we need to find, and then perform the calculation for, just the legs
  // in that path.
  vector <Path> paths;
  for (unsigned int j = 0; j < ScatteringPaths.size (); j ++) {
    if (ScatteringPaths[j].Index == ScatteringPaths[PathNo].Index) {
      paths.push_back (ScatteringPaths[j]);
    }
  }

  // Now obtain dChi for the selected scattering path and filter the result
  vector <Coordinate> dChij = doDChiCalculation (paths, getExptKData ());
  //for (unsigned int i = 0; i < dChij.size (); i ++) {
    //cout << PathNo << "  " << dChij[i].x << ", " << dChij[i].y << endl;
  //}
  FFData dChijData;
  fourierFilter (fkMin, fkMax, fdK, frMin, frMax, fdR, dChij, dChijData);
  return dChijData;
}


int ThermalSpectrum::copySpectrum (ThermalSpectrum Input) {
  DataPoints = Input.getDataPoints ();
  AllDataPoints = DataPoints;
  RawData = DataPoints;
  SigmaSqr = Input.getSigmaSqr ();
  return DX_NO_ERROR;
}
