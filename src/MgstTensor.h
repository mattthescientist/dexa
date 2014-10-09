// MgstTensor.h
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
#ifndef MT_MGSTTENSOR_H
#define MT_MGSTTENSOR_H

#define MT_NO_ERROR 0
#define MT_ERROR    1
#define MT_OUTPUT_PRECISION     5

#define MT_TENSOR_CUBIC         1
#define MT_TENSOR_CYLINDRICAL   2
#define MT_TENSOR_TETRAGONAL_T1 3
#define MT_TENSOR_TETRAGONAL_T2 4
#define MT_TENSOR_ISOTROPIC     5
#define MT_DEF_SYMMETRY         MT_TENSOR_CUBIC

#ifndef PI
  #define PI (4.0 * atan (1.0))
#endif

#include <vector>
#include "DiffXasSpectrum.h"

using namespace::std;

class MgstTensor {

public:
  // Constructors (x2) and destructor
  MgstTensor ();
  MgstTensor (double LA0, double LG2, double LE2, int Sym) : 
    theLambdaA0 (LA0), theLambdaG2 (LG2), theLambdaE2 (LE2), theSymmetry (Sym) 
    { updateVoigtMatrix (); updateMgstTensor (VoigtMatrix); }
  virtual ~MgstTensor () {};

  // Get and Set routines for the tensor properties
  virtual int mgstCoefficients (double NewLA0, double NewLG1, double NewLE2);
  virtual int mgstCoefficients (double NewLG1, double NewLD2, double NewL1A0,
    double NewL2A0, double NewL1A2, double NewL2A2);

  // Get set routines for magnetostriction coefficients. See private section
  // below for which coefficients correspond to which crystal symmetries
  virtual int lambdaA0 (double NewLambdaA0);
  virtual int lambdaG2 (double NewLambdaG2);
  virtual int lambdaE2 (double NewLambdaE2);
  virtual int lambdaD2 (double NewLambdaD2);
  virtual int lambda1A0 (double NewLambda1A0);
  virtual int lambda2A0 (double NewLambda2A0);
  virtual int lambda1A2 (double NewLambda1A2);
  virtual int lambda2A2 (double NewLambda2A2);
  virtual int lambda1G2 (double NewLambda1G2);
  virtual int lambda2G2 (double NewLambda2G2);
  virtual int lambda3G2 (double NewLambda3G2);
  virtual int lambda4G2 (double NewLambda4G2);  
  virtual int tensorElement (int Element, double NewValue);
  virtual int voigtElement (int Element, double NewValue);
  virtual int halfVoigtElement (int Element, double NewValue);
  
  virtual double lambdaA0 () { return theLambdaA0; }
  virtual double lambdaG2 () { return theLambdaG2; }
  virtual double lambdaE2 () { return theLambdaE2; }
  virtual double lambdaD2 () { return theLambdaD2; }
  virtual double lambda1A0 () { return theLambda1A0; }
  virtual double lambda2A0 () { return theLambda2A0; }
  virtual double lambda1A2 () { return theLambda1A2; }
  virtual double lambda2A2 () { return theLambda2A2; }
  virtual double lambda1G2 () { return theLambda1G2; }
  virtual double lambda2G2 () { return theLambda2G2; }
  virtual double lambda3G2 () { return theLambda3G2; }
  virtual double lambda4G2 () { return theLambda4G2; }
  virtual double tensorElement (int i) { 
    return Tensor 
      [int(double(i)/27.0)%3]
      [int(double(i)/9.0)%3]
      [int(double(i)/3.0)%3]
      [i%3]; 
  }
  virtual double voigtElement (int i);
  virtual double halfVoigtElement (int i);
  
  virtual int symmetry (int NewSymmetry);
  virtual int symmetry () { return theSymmetry; }

  vector < vector < double > > getVoigtMatrix ();
  vector < vector < vector < vector < double > > > > getMgstTensor ();

  int rotate (double x, double y, double z);
  virtual int prefOrientation (double x, double y, double z);

  // Some output routines for debugging purposes
  int saveVoigt (char *Filename);
  int saveTensor (char *Filename);
  int rotateDebug (double x, double y, double z);
  
private:
  int updateVoigtMatrix ();
  int updateMgstTensor (double theVoigtMatrix [6][6]);
  vector < vector <double> > bondMatrix(vector < vector <double> > DirCosines);
  int voigtPrefOrientation ();
  
  // Cubic coefficients (and isotropic if E2 is excluded)
  double theLambdaA0;
  double theLambdaG2;
  double theLambdaE2;
  
  // Additional coefficients for Cylindrical and Tetragonal T1
  double theLambda1A0, theLambda2A0;
  double theLambda1A2, theLambda2A2;
  double theLambdaD2;
  
  // Coefficients for Tetragonal T2 which are added onto T1 coefficients above
  // minus lambdaD2
  double theLambda1G2, theLambda2G2;
  double theLambda3G2, theLambda4G2;

  int theSymmetry;
    
protected:
  vector < vector <double> > eulerMatrix (double x, double y, double z);

  Coordinate PrefOrientation; // Mean orientation of crytals wrt beamline frame
  double VoigtMatrix [6][6];  // Reduced tensor notation
  double Tensor [3][3][3][3]; // 4th-rank Magnetostriction tensor
  double VoigtParams [6][6];  // Fit parameters for the Voigt Matrix;
};

#endif // MT_MGSTTENSOR_H
