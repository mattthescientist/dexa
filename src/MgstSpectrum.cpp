// MgstSpectrum.cpp
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
#include "MgstSpectrum.h"

#include <cmath>

// This class inherits and so completes the DiffXasSpectrum abstract base class
// by providing 'dChi' to calculate DiffXAS signals for MAGNETOSTRICTION.
//
// When a call is made to 'dChi' it first checks the 'SpectrumReady' flag to see
// if the information required to calculate the DiffXAS signal is present. If it
// is not, 'contractPaths' is called to calculate and store the atomic strain
// 'dS' for for each individual scattering path.
//
// In order for 'contractPaths' to do this correctly, the user must have set the
// following class properties before calling 'dChi'
//
//   1. 'polarisation'    : to specify the x-ray polarisation direction
//   2. 'prefOrientation' : to set the preferential orientation of crystalites
//   3. 'magnetisation1'  : to set the first sample magnetisation vector
//   4. 'magnetisation2'  : to set the second sample magnetisation vector
//   5. 'lambdaA0'        : to set the Lambda A0 magnetostriction coefficient
//   6. 'lambdaG2'        : to set the Lambda G2 magnetostriction coefficient
//   7. 'lambdaE2'        : to set the Lambda E2 magnetostriction coefficient
//   (5 - 7 assume cubic symmetry. If this differs set the appropriate coeffs.)
//
// If these are not set, 'contractPaths' will just work with the default values.
// Given the default magnetostriction coefficients are zero, this will not
// produce very interesting results!
//
// When any of the above routines are called, 'SpectrumReady' is set to false,
// causing all the magnetostrictive strains to be recalculated on the next call
// to 'dChi'. Equally, once 'contractPaths' has executed successfully,
// 'SpectrumReady' is set to true, indicating that the strains 'dS' are ready to
// be read out.
//
// 'contractPaths' evaluates the mean magnetostrictive strain for each 
// scattering path based on the assumption that the polycrystals in the sample 
// are preferentially oriented in some known manner along the line of the beam 
// and randomly oriented in the plane perpendicular to the beam. It thus works
// out the mean strain resulting from a configurational average of the 
// magnetostriction about the line of the beam. The 'numSteps' property may be 
// used to change the number of steps at which to evaluate the strains about 
// this axis.
//
// The strain is calculated for each leg explicitly by using 'contractLegs', the
// results from each leg being summed to obtain the total path strain. For each
// leg, 'contractToTensor' is called to calculate
//
//   dS_{leg} =  T_{ijkl} * l_{i} * l_{j} * M_{k} * M_{l}
//
// where T is the magnetostriction tensor, l the leg vector, and M the
// magnetisation vector. Needless to say that 'contractToTensor' is called twice
// for each leg so that dS may be obtained for each magnetisation direction.
//
// The final job for 'contractLegs' is to identify the first and last leg in any
// scattering path and additionally contract them onto the x-ray polarisation
// vector. This is what actually gives rise to the dichroic signal.
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
void MgstSpectrum::calculateDChi () {

  vector <Coordinate> temp;
  vector <Coordinate> dChi = doDChiCalculation (dS, temp);

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
  //for (unsigned int i = 0; i < TheoryFFData.QData.size (); i ++) {
    //cout << "Q " << i << " " << TheoryFFData.QData[i].x << ", " << TheoryFFData.QData[i].y << endl;
  //}
  generateDChiCheb (dChiSampled);  
}


//------------------------------------------------------------------------------
// getDChiPath (int) : returns dChi_j(k). That is, the contribution to the total
// dChi signal by path j only. This can be used for diagnostic/debug purposes to
// see the effect of individual paths on the overall dChi. As with 
// calculateDChi() above, the actual calculation is done in the sub-function
// doDChiCalculation() below. The difference is that only the scattering legs of
// path j are passed to it rather than all the legs.
//
FFData MgstSpectrum::getDChiPath (int PathNo) {

  // First make sure that all the strains have been calculated. If
  // they haven't, or need updating, call 'contractPaths' to do that.
  if (!SpectrumReady) {
    calculateMagnetostriction ();
    SpectrumReady = true;
  }
  
  // Since we are only returning the dChi component for a single scattering 
  // path, we need to find, and then perform the calculation for, just the legs
  // in that path.
  vector <DeltaS> thePaths;
  for (unsigned int j = 0; j < dS.size (); j ++) {
    if (dS[j].Chip -> Index == ScatteringPaths[PathNo].Index) {
      thePaths.push_back (dS[j]);
      //cout << "DchiPath " << j << " " << PathNo << " " << dS[j].Chip -> Index << endl;
    }
  }

  // Now obtain dChi for the selected scattering path and filter the result
  vector <Coordinate> dChij = doDChiCalculation (thePaths, getExptKData ());
  //for (unsigned int i = 0; i < dChij.size (); i ++) {
    //cout << PathNo << "  " << dChij[i].x << ", " << dChij[i].y << endl;
  //}
  FFData dChijData;
  fourierFilter (fkMin, fkMax, fdK, frMin, frMax, fdR, dChij, dChijData);
  return dChijData;
}


//------------------------------------------------------------------------------
// doDChiCalculation (...) : This function is a wrapper for doDchiCalculation1
// below. It allows dChi to be calculated at either regular intervals of
// EXAFS_K_STEP between EXAFS_K_MIN and EXAFS_K_MAX, or at specific points. If
// XCoords is blank at the input, it will be the former, if not, it will be the
// latter, where the points will be taken from the XCoords.x property.
//
vector <Coordinate> MgstSpectrum::doDChiCalculation (vector <DeltaS> &dSIn, vector <Coordinate> XCoords) {
  vector <Coordinate> dChi;

  double k, kIn;
  if (XCoords.size () == 0) {
    for (kIn = EXAFS_MIN_K; k < EXAFS_MAX_K; kIn += EXAFS_K_STEP) {    
      // Only proceed with the calculation if k > 0.0 after any edge shifts
      if (((kIn * kIn * ALPHA) - dE()) / ALPHA - pow(dK(),2) > 0.0) {
        k = pow(((kIn * kIn * ALPHA) - dE()) / ALPHA, 0.5) - dK();
        
        // Calculate dChi here. Pass in kIn, rather than k, as that must be
        // stored in the x ordinate of the output data points.
        dChi.push_back (doDChiCalculation1 (dSIn, kIn));
      }
    }
  } else {
    for (unsigned int j = 0; j < XCoords.size (); j ++) {
      k = XCoords[j].x;
      dChi.push_back (doDChiCalculation1 (dSIn, k));
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
Coordinate MgstSpectrum::doDChiCalculation1 (vector <DeltaS> &dSIn, double kIn){

  double k;
  Coordinate NewDChi;
  vector <Coordinate> dChi;   // dChi(k) evaluated at each k in the feffnnnn's

  // Calculate the DiffXAS at each k by summing contributions over each path j
  // Take the edge shift parameters into account when calculating the value of
  // k to use for the calculation. Take special care to check that any dE does
  // not result in us calculating dChi for k <= 0.0
  NewDChi.y = 0.0;
  if (((kIn * kIn * 3.81) - dE()) / 3.81 - pow(dK(),2) > 0.0) {
    k = pow(((kIn * kIn * ALPHA) - dE()) / ALPHA, 0.5) - dK();
    for (unsigned int j = 0; j < dSIn.size (); j ++) {

      //----------------------------------------------------------------------
      // Standard EXAFS components. These are the data read from the feffnnnn
      // files. See DiffXasSpectrum::readFeffFile(...) for more information.
      //
      double A    = gsl_cheb_eval (dSIn[j].Chip -> ACheb, k)*dSIn[j].Chip ->S02;
      double P    = gsl_cheb_eval (dSIn[j].Chip -> PCheb, k);
      double Phi  = gsl_cheb_eval (dSIn[j].Chip -> PhiCheb, k);
      double Lam  = gsl_cheb_eval (dSIn[j].Chip -> LamCheb, k);
      double Sig2 = dSIn[j].Chip -> Sig2;
      double S    = dSIn[j].Chip -> S;

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
      // dSj is the change in scattering path length due to magnetostriction.
      // Note that dS[j].Ave is the NORMALISED change in path length. i.e. it
      // is the strain. It is therefore multiplied by the path length to get
      // the actual atomic displacement.
      //
      // Cj is a COMMON amplitude factor, which contains all the components of
      // the EXAFS function that do not change between measurements the
      // magnetisation measurements at M' and M''.
      //
      // Dj is the amplitude factor that DEPENDS upon the magnetisation, and
      // so changes between the two measurements at M' and M''.
      //
      // Qj is the phase component of EXAFS signal. Like Dj, it is dependent
      // upon the magnetisation, and so changes between the two measurements.
      //
      // Each variable is assigned the parameters for M', the first 
      // magnetisation state, upon declaration.
      //
      double Rj  = S * 0.5;
      double dRj = dSIn[j].Chip -> dR;
      double dSj = dSIn[j].Ave * (S + 2.0 * dRj);
      double Cj  = A * pow (Rj, 2) * exp (-2.0 * P * P * Sig2 
        + 2.0 * Sig2 / (Lam * Lam) - 2.0 * dRj / Lam + 4.0 * Sig2 / (Rj * Lam));
      double Dj  = pow (dRj + Rj, -2.0);
      double Qj  = (2.0 * k * Rj) + Phi + (2.0 * P * dRj) - (4 * P * Sig2 / Rj)
        - (4.0 * P * Sig2 / Lam);


      // These commented out lines can be inserted to produce dChi based upon
      // the Taylor approximation formalism.
      // double dDj = 1.0 / (2.0 * (dRj + Rj));
      // double dQj = 2.0 * P;
      // NewDChi.y += Cj * ( Dj * dQj * cos (Qj) + dDj * sin (Qj)) * dSj / 2.0;
      // NewDChi.y += Cj * Dj * cos (Qj) * k * dSj;
      // NewDChi.y += A * pow (Rj, 2) * exp (-2.0 * k * k * Sig2) * Dj * cos (2.0 * k * Rj + Phi) * k * dSj;

      // Now change the magnetisation to M'' and recalculate the dependent
      // components of the DiffEXAFS function (i.e. Dj and Qj)
      dRj += dSj / 2;
      double Dj2 = pow (dRj + Rj, -2.0);
      double Qj2 = (2.0 * k * Rj) + Phi + (2.0 * P * dRj) - (4 * P * Sig2 / Rj)
        - (4.0 * P * Sig2 / Lam);
      // Finally, take the difference between the two measurements to get the
      // overall DiffEXAFS signal.
      NewDChi.y += Cj * ( (Dj * sin (Qj)) - (Dj2 * sin (Qj2)) );
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
void MgstSpectrum::generateDChiCheb (vector <Coordinate> &Data) {
  gsl_cheb_free (dChiCheb);
  dChiCheb = gsl_cheb_alloc (CHEB_ORDER);
  gsl_function F;
  F.params = &Data;
  F.function = &MgstSpectrum::dChiChebFunction;
  gsl_cheb_init (dChiCheb, &F, Data[0].x, Data[Data.size() - 1].x);
}


//------------------------------------------------------------------------------
// dChiChebFunction (double, void *) : Returns dChi(k) for any arbirary k by
// using linear interpolation. This function is used by GSL when generating a
// Chebyshev approximation of dChi. Thereafter, the Chebyshev is used to obtain 
// dChi(k) rather than this function.
//
double MgstSpectrum::dChiChebFunction (double k, void *p) {
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
double MgstSpectrum::dChi (double k) {

  // First make sure that all the strains have been calculated. If
  // they haven't, or need updating, call 'contractPaths' to do that.
  if (!SpectrumReady) {
    calculateMagnetostriction ();
    SpectrumReady = true;
  }
  // Then obtain dChi(k) from the dChiCheb Chebyshev
  if (((k * k * 3.81) - dE()) / 3.81 - pow(dK(),2) <= 0.0) return 0.0;
  return gsl_cheb_eval (dChiCheb, k);
}


//------------------------------------------------------------------------------
// calculateMagnetostriction () : If a call is made to either of the two chi() 
// or dChi() functions and the SpectrumReady flag is false, this function is 
// called upon to calculate the magnetostriction signal components. This
// particular function is responsible for orientating the crystal along the
// previously specified preferential orientation directions and then calling 
// contractPaths() to obtain the magnetostrictive strains. 
//
// If the contribution from given preferential directions sum to less than unity
// it is assumed that the remainder of crystallites are randomly orientated.
// contractPaths() is then called for each of NumAveSteps orientations between 0
// and PI radians about the beamline frame's x-axis. When this is combined with
// the 0 to 2PI averaging around the beam propagation vector in contractPaths(),
// the resulting average is that taken over a whole sphere.
//
void MgstSpectrum::calculateMagnetostriction () {
  Coordinate theAtom;
  vector<DeltaS> dSAveraged;
  double DegreeRandom = 1.0;

  // Perform the configurational average of the scattering paths around the
  // specified preferential directions. Weight each direction according to the
  // fraction of crystallites orientated along it (i.e the degree of orientation
  // in that direction).
  for (unsigned int i = 0; i < PrefOrientations.size (); i ++) {
    for (unsigned int j = 0; j < Atoms.size (); j ++) {
      theAtom.x = OriginalAtoms[j].x; 
      theAtom.y = OriginalAtoms[j].y; 
      theAtom.z = OriginalAtoms[j].z;
      rotateVector (theAtom, 
        PrefOrientations[i].x, PrefOrientations[i].y, PrefOrientations[i].z);
      Atoms[j].x = theAtom.x; Atoms[j].y = theAtom.y; Atoms[j].z = theAtom.z; 
    }
    MgstTensor::prefOrientation 
      (PrefOrientations[i].x, PrefOrientations[i].y, PrefOrientations[i].z);
    contractPaths ();
    dSAveraged.resize (dS.size ());
    for (unsigned int j = 0; j < dS.size (); j ++) {
      dSAveraged[j].AveMag1 += dS[j].AveMag1 * PrefOrientations[i].Degree;
      dSAveraged[j].AveMag2 += dS[j].AveMag2 * PrefOrientations[i].Degree;
      dSAveraged[j].Ave += dS[j].Ave * PrefOrientations[i].Degree;
    }
    DegreeRandom -= PrefOrientations[i].Degree;
  }
  
  // Now assess the degree of random orientation in the crystallites and average
  // the scattering path contributions over all angles.
  if (DegreeRandom > 0.0) {
    for (unsigned int Step = 0; Step < NumAveSteps; Step ++) {
      for (unsigned int j = 0; j < Atoms.size (); j ++) {
        theAtom.x = OriginalAtoms[j].x; 
        theAtom.y = OriginalAtoms[j].y; 
        theAtom.z = OriginalAtoms[j].z;
        rotateVector 
          (theAtom, double(Step) * 2.0 * PI / double (NumAveSteps), 0.0, 0.0);
        Atoms[j].x = theAtom.x; Atoms[j].y = theAtom.y; Atoms[j].z = theAtom.z; 
      }
      MgstTensor::prefOrientation 
        (double(Step) * 2.0 * PI / double (NumAveSteps), 0.0, 0.0);
      contractPaths ();
      dSAveraged.resize (dS.size ());
      for (unsigned int j = 0; j < dS.size (); j ++) {
        dSAveraged[j].AveMag1 += dS[j].AveMag1 * DegreeRandom / NumAveSteps;
        dSAveraged[j].AveMag2 += dS[j].AveMag2 * DegreeRandom / NumAveSteps;
        dSAveraged[j].Ave += dS[j].Ave * DegreeRandom / NumAveSteps;
      }
    }
  }
  
  // Finally, replace the dS vector with the dSAveraged vector created here.
  for (unsigned int i = 0; i < dS.size (); i ++) {
    dS[i].AveMag1 = dSAveraged[i].AveMag1;
    dS[i].AveMag2 = dSAveraged[i].AveMag2;
    dS[i].Ave = dSAveraged[i].Ave;
  }
  
  calculateDChi ();
}


//------------------------------------------------------------------------------
// contractPaths () : Goes through all the loaded scattering paths, one shell at
// a time, and calculates the mean dS due to Magnetostriction, by performing an
// orientational average of sample crystalites about the axis of the beam 
// propagation direction.
//
void MgstSpectrum::contractPaths () {
  vector <DeltaS*> *newdS;
  DeltaS currentdS;
  int CurrentPath;
  
  // Handle the first orientation explicitly so as to generate the dS vector
  dS.erase (dS.begin(), dS.end());
  rotate (0.0, 0.0, 0.0);
  for (unsigned int theShell = 0; theShell < Shells.size (); theShell ++) {
    if (getPath (Shells[theShell].ChipIndex) != NULL) { 
      newdS = new vector <DeltaS*>;
      Coordinate Blank;
      contractLegs (Paths[theShell], 
        0.0, 0.0, theShell, 0.0, 0.0, 0.0, 0.0, Blank, newdS);
      for (unsigned int i = 0; i < newdS -> size (); i ++) {
        currentdS.AveMag1 = newdS -> at (i) -> AveMag1 / NumAveSteps;
        currentdS.AveMag2 = newdS -> at (i) -> AveMag2 / NumAveSteps;
        currentdS.Ave = newdS -> at (i) -> Ave / NumAveSteps;
        currentdS.Chip = newdS -> at (i) -> Chip;
        dS.push_back (currentdS);
      }
      for (unsigned int i = 0; i < newdS -> size (); i ++) { delete newdS -> at(i); }
      delete newdS;
    } //else { cout << "Null chipindex at " << theShell << endl; }
  }
  
  // Now repeat the contraction process for each tensor orientation and update
  // the average dS values for each path.
  for (unsigned int theStep = 1; theStep < NumAveSteps; theStep ++) {
    double Theta = double(theStep) * 2.0 * PI / double(NumAveSteps);
    rotate (0.0, 0.0, Theta);
    CurrentPath = 0;
    for (unsigned int theShell = 0; theShell < Shells.size (); theShell ++) {
      if (getPath (Shells[theShell].ChipIndex) != NULL) {
        newdS = new vector <DeltaS*>;
        Coordinate Blank;
        contractLegs(Paths[theShell], 
          0.0, 0.0, theShell, 0.0, 0.0, Theta, 0.0, Blank, newdS);
        for (unsigned int i = 0; i < newdS -> size (); i ++) {
          dS[i + CurrentPath].AveMag1 += newdS -> at (i) -> AveMag1/NumAveSteps;
          dS[i + CurrentPath].AveMag2 += newdS -> at (i) -> AveMag2/NumAveSteps;
          dS[i + CurrentPath].Ave += newdS -> at (i) -> Ave / NumAveSteps;
        }
        CurrentPath += newdS -> size ();
        for (unsigned int i = 0; i < newdS -> size (); i ++) { delete newdS -> at(i); }
        delete newdS;
      }
    }
  }
  
/*  for (unsigned int i = 0; i < ScatteringPaths.size (); i ++) ScatteringPaths[i].dSj.clear ();
  for (unsigned int i = 0; i < dS.size (); i ++) {
//    dS[i].Ave += dS[i].Chip -> dS;
    dS[i].Chip -> dSj.push_back (dS[i].Ave);
  }*/
}


//------------------------------------------------------------------------------
// rotateVector (Coordinate *, double, double, double) : Generates the Euler
// rotation matrix from the three parameters x, y, z - which specify rotations
// about those respective axes in RADIANS - and then rotates theVector
// accordingly. As always, the coordinates are with respect to the beamline
// frame.
//
int MgstSpectrum::rotateVector (
  Coordinate &theVector, double x, double y, double z) {
  
  vector < vector <double> > EulerMatrix = eulerMatrix (x, y, z);
  double theOldVector [3];
  double theNewVector [3];
  
  theOldVector [0] = theVector.x;
  theOldVector [1] = theVector.y;
  theOldVector [2] = theVector.z;
  for (int i = 0; i < 3; i ++) {
    theNewVector [i] = 0.0;
    for (int j = 0; j < 3; j ++) {
      theNewVector [i] += EulerMatrix [i][j] * theOldVector [j];
    }
  }
  theVector.x = theNewVector [0];
  theVector.y = theNewVector [1];
  theVector.z = theNewVector [2];
  return MS_NO_ERROR;
}
  
  
//------------------------------------------------------------------------------
// prefOrientation (double, double, double) : Rotates the crystal lattice to
// preferential orientation with respect to the beamline frame. This doesn't 
// just rotate the atomic positions, but also the Magnetostriction tensor. Both
// these rotations are absolute rather than relative, and are with respect to
// the beamline frame.
//
int MgstSpectrum::addPrefOrientation (double x,double y,double z,double Degree){
  if (Degree > 0.0) {
    double DegreeSoFar = 0.0;
    for (unsigned int i = 0; i < PrefOrientations.size (); i ++) {
      DegreeSoFar += PrefOrientations[i].Degree;
    }
    if (DegreeSoFar + Degree > 1.0) {
      cout << "Cannot add preferential orientation (" << x << ", " << y << ", "
        << z << ") with a degree of " << Degree * 100.0 << "%.\n The total "
        << "degree of preferential orientation over all specified directions "
        << "must sum to less than or equal to 1.0." << endl;
      return MS_ERROR;
    }
    PrefOrient NewPrefOrientation;
    NewPrefOrientation.x = x;
    NewPrefOrientation.y = y;
    NewPrefOrientation.z = z;
    NewPrefOrientation.Degree = Degree;
    PrefOrientations.push_back (NewPrefOrientation);
    SpectrumReady = false;
    return MS_NO_ERROR;
  }
  cout << "Invalid preferential orientation. Degree must be greater than zero." 
    << endl;
  return MS_ERROR;
} 

int MgstSpectrum::prefX (int Orientation, double NewX) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    PrefOrientations[Orientation].x = NewX;
    return MS_NO_ERROR;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return MS_ERROR; 
  }
}

int MgstSpectrum::prefY (int Orientation, double NewY) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    PrefOrientations[Orientation].y = NewY;
    return MS_NO_ERROR;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return MS_ERROR; 
  }
}
int MgstSpectrum::prefZ (int Orientation, double NewZ) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    PrefOrientations[Orientation].z = NewZ;
    return MS_NO_ERROR;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return MS_ERROR; 
  }
}
int MgstSpectrum::prefDeg (int Orientation, double NewDeg) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    if (NewDeg < 0.0) PrefOrientations[Orientation].Degree = 0.0;
    else if (NewDeg > 1.0) PrefOrientations[Orientation].Degree = 1.0;
    else PrefOrientations[Orientation].Degree = NewDeg;
    return MS_NO_ERROR;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return MS_ERROR; 
  }
}

double MgstSpectrum::prefX (int Orientation) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    return PrefOrientations[Orientation].x;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return 0.0; 
  }
}
double MgstSpectrum::prefY (int Orientation) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    return PrefOrientations[Orientation].y;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return 0.0; 
  }
}
double MgstSpectrum::prefZ (int Orientation) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    return PrefOrientations[Orientation].z;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return 0.0; 
  }
}
double MgstSpectrum::prefDeg (int Orientation) {
  if (Orientation >= 0 && Orientation < int(PrefOrientations.size ())) {
    return PrefOrientations[Orientation].Degree;
  } else {
    cout << "Invalid orientation index: " << Orientation << endl;
    return 0.0; 
  }
}
  
//------------------------------------------------------------------------------
// contractLegs (vector <Leg*>, double, double, int, vector <DeltaS>) :
// Contracts each leg of the paths in arg1 onto the Magnetostriction tensor. The
// additional contractions due to x-ray polarisation are calculated from the 
// first and last legs in each path and are responsible for the dichroic signal.
//
// In order for this routine to contract every leg in the linked list of path
// legs, it is called RECURSIVELY; using the Next pointer to move from one leg
// to the following. When the Next pointer is NULL, the end of the path has
// been reached and the total dS for the path pushed into a vector of all paths.
//
void MgstSpectrum::contractLegs (vector <Leg*> *PathLegs, double dSMag1In, 
  double dSMag2In, int theShell, double x, double y, 
  double z, double Pol1, Coordinate PrevLeg, vector <DeltaS*> *alldS) {
  double dSMag1, dSMag2;
  double VectorLength;
  
  for (unsigned int i = 0; i < PathLegs -> size (); i ++) {
    if (PathLegs -> at (i) -> Next != NULL) {
    
      // First generate the unit vector for the leg and get its length. Also
      // rotate the leg to the orientation specified by the x, y, and z input
      // arguments.
      Coordinate theLeg;
      theLeg.x = Atoms[PathLegs -> at (i) -> Next -> at (0) -> Atom].x
        - Atoms [PathLegs -> at (i) -> Atom].x;
      theLeg.y = Atoms[PathLegs -> at (i) -> Next -> at (0) -> Atom].y
        - Atoms [PathLegs -> at (i) -> Atom].y;
      theLeg.z = Atoms[PathLegs -> at (i) -> Next -> at (0) -> Atom].z
        - Atoms [PathLegs -> at (i) -> Atom].z;
      VectorLength = vectorLength (theLeg.x, theLeg.y, theLeg.z);
      theLeg.x /= VectorLength;
      theLeg.y /= VectorLength;
      theLeg.z /= VectorLength;
      
      rotateVector (theLeg, x, y, z);
      //cout << theLeg.x << "," << theLeg.y << "," << theLeg.z << endl;
      
      // Add the dS for this leg to the total for the path so far. This is done
      // twice - once for each magnetisation direction.
      dSMag1 = dSMag1In + /*VectorLength **/contractToTensor(theLeg,Magnetisation1);
      dSMag2 = dSMag2In + /*VectorLength **/contractToTensor(theLeg,Magnetisation2);
      
      // If the current leg is the first in the path, a polarisation factor
      // must also be calculated, which is used later.
      if (dSMag1In == 0.0 && dSMag2In == 0.0) {
        Pol1 = abs (theLeg.x * Polarisation.x 
          + theLeg.y * Polarisation.y + theLeg.z * Polarisation.z);
      }
      
      // Move on to the next leg in the path, sending on dSMag1 and dSMag2 - the
      // cumulative dS's for every leg consider so far.
      contractLegs (PathLegs -> at(i) -> Next, 
        dSMag1, dSMag2, theShell, x, y, z, Pol1, theLeg, alldS);
    } else {

      // The end of the path has been reached. Since the current leg is the 
      // last, generate the second of the two polarisation factors.
      double Pol2 = abs (PrevLeg.x * Polarisation.x 
        + PrevLeg.y * Polarisation.y + PrevLeg.z * Polarisation.z);
      
      // Finally, complete the calculation for the two dS's and take their also
      // store the average of the two. Push this information onto thedS vector,
      // which is now returned to the 'contractPaths' routine above.
      DeltaS *thedS = new DeltaS;
      thedS -> AveMag1 = dSMag1In * Pol1 * Pol2 /*/ Shells[theShell].Degeneracy*/;
      thedS -> AveMag2 = dSMag2In * Pol1 * Pol2 /*/ Shells[theShell].Degeneracy*/;
      thedS -> Ave = thedS -> AveMag2 - thedS -> AveMag1;
      thedS -> Chip = getPath(Shells[theShell].ChipIndex);
      alldS -> push_back (thedS);
      //cout << dSMag1In << " " << dSMag2In << " " << Pol1 << " " << Pol2 << "      " << thedS -> AveMag1 << "  " << thedS -> AveMag2 << "  " << thedS->Ave << endl;
    }
  }
}


//------------------------------------------------------------------------------
// contractToTensor (Coordinate, double [3][3]) : Contracts the input scattering
// leg and magnetisation vectors onto the magnetostriction tensor (see 
// MgstTensor.cpp for more on the tensor itself).
//     
double MgstSpectrum::contractToTensor (
  Coordinate theLegIn, Coordinate Magnetisation1) {
  double theLeg [3];
  theLeg [0] = theLegIn.x; theLeg [1] = theLegIn.y; theLeg [2] = theLegIn.z;

  double Mag[3];
  Mag [0] = Magnetisation1.x;
  Mag [1] = Magnetisation1.y;
  Mag [2] = Magnetisation1.z;


  double Sum = 0.0;
  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 3; j ++) {
      for (int k = 0; k < 3; k ++) {
        for (int l = 0; l < 3; l ++) {
          Sum += Tensor [i][j][k][l] * theLeg[i] * theLeg[j] * Mag[k] * Mag[l];
          //cout << i << "," << j << "," << k << "," << l << ": " << Tensor[i][j][k][l] << "  " << theLeg[i] << "  " << theLeg[j] << "  " << Mag[k] << "  " << Mag[l] << "  " << Sum << endl;
        }
      }
    }
  }
  return Sum;
}


int MgstSpectrum::copySpectrum (MgstSpectrum Input) {
  DataPoints = Input.getDataPoints ();
  AllDataPoints = DataPoints;
  RawData = DataPoints;
  SigmaSqr = Input.getSigmaSqr ();
  return DX_NO_ERROR;
}
