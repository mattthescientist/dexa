// MgstTensor.cpp
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
#include "MgstTensor.h"

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>

using namespace::std;

// The MgstTensor class provides all the necessary functionality for creating
// and manipulating the 4th-rank Magnetostriction tensor. At the moment, this is
// restricted to crystals with CUBIC SYMMETRY ONLY. The rules for generating the
// tensor based on other symmetries must be added to this code if so desired.
//
// Tensor Generation:
// -----------------
// To use this class, the user must either set the 'mgstCoefficients' property -
// passing the values of the lambda-alpha-zero, lambda-gamma-two, and 
// lambda-epsilon-two coefficients as arguments - or use the properties 
// 'lambdaA0', 'lambdaG2', and 'lambdaE2' to set each coefficient independently.
// Each time one of these properties is changed, the tensor is reconstructed via
// contracted Voigt notation.
//
// This is done in two steps. First, 'updateVoigtMatrix' is called to construct
// the 6x6 element Voigt matrix based on the input coefficient values. The rules
// for generating this matrix may be found in "The Physical Properties of
// Crystals" by J.F. Nye. Once this matrix is complete, 'updateMgstTensor' is
// called to expand the 6x6 Voigt matrix into the full 3x3x3x3 magnetostriction
// tensor.
//
// At any time, the user may call 'getVoigtMatrix' or 'getMgstTensor' to obtain
// a copy of either the Voigt matrix or magnetostriction tensor respectively.
// Both are returned by value in the form of a std::vector array, so the user 
// may subsequently alter or even delete them without affecting the originals
// stored in this class.
//
// In reality though, the user should not need to use either the Voigt matrix or
// magnetostriction tensor directly since this class is inherited by 
// MgstSpectrum. Consequently, beyond the need to set the magnetostriction 
// coefficients in the first place, the user will just interact with 
// MgstSpectrum to obtain the DiffXAS information. 
//
// Tensor Rotation:
// ---------------
// It is important to note that none of the rotation routines here operate on
// the magnetostriction tensor directly. Instead, for increased speed and
// computational efficiency, the rotation is achieved by using Bond transform
// matricies to rotate the contracted Voigt matrix. Once rotated, the Voigt
// matrix is then expanded to the full tensor via 'updateMgstTensor'
//
// The tensor may be rotated via either the 'rotate' method, or the 
// 'prefOrientation' property. In both cases, three arguments must be provided // for rotation IN RADIANS about the beamline's x, y, and z axes 
// respectively.

// The two routines are essentially the same but differ in one key respect:
// 'rotate' does not permanently store the rotated Voigt matrix whereas
// 'prefOrientation' OVERWRITES the original Voigt matrix. Consequently,
// 'prefOrientation' is used to set a reference orientation (with respect to the
// beamline frame) with 'rotate' then being to further rotate the crystal with
// respect to this reference. Consecutive calls to 'rotate' are not cumulative -
// each rotation is always with respect to the 'prefOrientation' reference.
//
// Note that the preferential orientation IS LOST when any of the 
// magnetostriction coefficients are changed. The user MUST therefore call 
// 'prefOrientation' after changing any coefficients.
//
//------------------------------------------------------------------------------
// Default constructor : Sets the crystal symmetry and magnetostriction
// coefficients to their default values as declared in MgstTensor.h.
// 
MgstTensor::MgstTensor () {
  theLambdaA0 = 0.0;
  theLambdaG2 = 0.0;
  theLambdaE2 = 0.0;
  theLambdaD2 = 0.0;
  theLambda1A0 = 0.0;
  theLambda2A0 = 0.0;
  theLambda1A2 = 0.0;
  theLambda2A2 = 0.0;
  theSymmetry = MT_DEF_SYMMETRY;
  PrefOrientation.x = 0.0; 
  PrefOrientation.y = 0.0; 
  PrefOrientation.z = 0.0; 
  for (int i = 0; i < 6; i ++) {
    for (int j = 0; j < 6; j ++) {
      VoigtParams [i][j] = 0.0;
      VoigtMatrix [i][j] = 0.0;
    }
  }
}


//------------------------------------------------------------------------------
// mgstCoefficients (double, double, double) : Sets all three magnetostriction
// coefficients at once. The arguments are for the lambda-alpha-zero, 
// lambda-gamma-two, and lambda-epsilon-two coefficients respectively. The 
// tensor is then regenerated.
//
int MgstTensor::mgstCoefficients (double NewLA0, double NewLG2, double NewLE2) {
  theLambdaA0 = NewLA0;
  theLambdaG2 = NewLG2;
  theLambdaE2 = NewLE2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

int MgstTensor::mgstCoefficients (double NewLG2, double NewLD2, double NewL1A0, 
  double NewL2A0, double NewL1A2, double NewL2A2) {
  theLambdaA0 = NewLG2;
  theLambdaG2 = NewLD2;
  theLambda1A0 = NewL1A0;
  theLambda2A0 = NewL2A0;
  theLambda1A2 = NewL1A2;
  theLambda2A2 = NewL2A2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambdaA0 (double) : Changes the lambda-alpha-zero coefficient and recreates
// the magnetostriction tensor accordingly.
//
int MgstTensor::lambdaA0 (double NewLA0) {
  theLambdaA0 = NewLA0;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambdaG2 (double) : Changes the lambda-gamma-two coefficient and recreates
// the magnetostriction tensor accordingly.
//
int MgstTensor::lambdaG2 (double NewLG2) {
  theLambdaG2 = NewLG2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambdaE2 (double) : Changes the lambda-epsilon-two coefficient and recreates
// the magnetostriction tensor accordingly.
//
int MgstTensor::lambdaE2 (double NewLE2) {
  theLambdaE2 = NewLE2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambdaD2 (double) : Changes the lambda-delta-two coefficient and recreates
// the magnetostriction tensor accordingly.
//
int MgstTensor::lambdaD2 (double NewLD2) {
  theLambdaD2 = NewLD2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambda1A0 (double) : Changes the lambda-one-alpha-zero coefficient and 
// recreates the magnetostriction tensor accordingly.
//
int MgstTensor::lambda1A0 (double NewL1A0) {
  theLambda1A0 = NewL1A0;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambda2A0 (double) : Changes the lambda-two-alpha-zero coefficient and 
// recreates the magnetostriction tensor accordingly.
//
int MgstTensor::lambda2A0 (double NewL2A0) {
  theLambda2A0 = NewL2A0;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambda1A2 (double) : Changes the lambda-one-alpha-two coefficient and 
// recreates the magnetostriction tensor accordingly.
//
int MgstTensor::lambda1A2 (double NewL1A2) {
  theLambda1A2 = NewL1A2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// lambda2A2 (double) : Changes the lambda-two-alpha-two coefficient and 
// recreates the magnetostriction tensor accordingly.
//
int MgstTensor::lambda2A2 (double NewL2A2) {
  theLambda2A2 = NewL2A2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

int MgstTensor::lambda1G2 (double NewL1G2) {
  theLambda1G2 = NewL1G2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

int MgstTensor::lambda2G2 (double NewL2G2) {
  theLambda2G2 = NewL2G2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

int MgstTensor::lambda3G2 (double NewL3G2) {
  theLambda3G2 = NewL3G2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

int MgstTensor::lambda4G2 (double NewL4G2) {
  theLambda4G2 = NewL4G2;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

int MgstTensor::tensorElement (int Element, double NewValue) {
  Tensor 
    [int(double(Element)/27.0)%3]
    [int(double(Element)/9.0)%3]
    [int(double(Element)/3.0)%3]
    [Element%3] = NewValue;
  return MT_NO_ERROR;
}

double MgstTensor::voigtElement (int i) {
  return VoigtMatrix [int(double(i)/6.0)%6][i%6];
}

int MgstTensor::voigtElement (int Element, double NewValue) {
  VoigtParams [int(double(Element)/6.0)%6][Element%6] = NewValue;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

double MgstTensor::halfVoigtElement (int Element) {
  int Row = int ((pow(8.0 * Element + 1, 0.5) - 1) / 2.0);
  int Col = (Element - (Row * (Row + 1)) / 2) % 6;
  return VoigtMatrix [Row][Col];
}

int MgstTensor::halfVoigtElement (int Element, double NewValue) {
  int Row = int ((pow(8.0 * Element + 1, 0.5) - 1) / 2.0);
  int Col = (Element - (Row * (Row + 1)) / 2) % 6;
  VoigtParams [Row][Col] = NewValue;
  VoigtParams [Col][Row] = NewValue;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}


//------------------------------------------------------------------------------
// symmetry (int) : Changes the crystal symmetry of the specimen, reforms the
// Voigt matrix, and recreates the magnetostriction tensor accordingly.
//
int MgstTensor::symmetry (int NewSymmetry) {
  theSymmetry = NewSymmetry;
  updateVoigtMatrix ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// prefOrientation (double, double, double) : Changes the preferential 
// orientation of the sample crystallites with respect to the beamline frame.
// All rotations are absolute and so are not cumulative over several calls.
//
int MgstTensor::prefOrientation (double x, double y, double z) { 
  PrefOrientation.x = x; 
  PrefOrientation.y = y; 
  PrefOrientation.z = z; 
  voigtPrefOrientation ();
  updateMgstTensor (VoigtMatrix);
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// getVoigtMatrix () : Returns the Voigt matrix in the form of a 2-D vector 
// of 6x6 elements. The matrix is returned by value so that any subsequent
// changes made to it do not affect the matrix stored in this class.
//
vector < vector <double> > MgstTensor::getVoigtMatrix () {
  vector < vector <double> > MatrixOut;
  vector <double> Row;
  for (int i = 0; i < 6; i ++) {
    for (int j = 0; j < 6; j ++) {
      Row.push_back (VoigtMatrix[i][j]);
    }
    MatrixOut.push_back (Row);
    Row.clear ();
  }
  return MatrixOut;
}

//------------------------------------------------------------------------------
// getMgstTensor () : Returns the magnetostriction tensor in the form of a 4-D 
// vector of 3x3x3x3 elements. The matrix is returned by value so that any 
// subsequent changes made to it do not affect the tensor stored in this class.
//
vector < vector < vector < vector <double> > > > MgstTensor::getMgstTensor () {
  vector < vector < vector < vector < double> > > > TensorOut;
  vector < vector < vector < double > > > RowA;
  vector < vector < double > > RowB;
  vector < double > RowC;
  
  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 3; j ++) {
      for (int k = 0; k < 3; k ++) {
        for (int l = 0; l < 3; l ++) {
          RowC.push_back (Tensor[i][j][k][l]);
        }
        RowB.push_back (RowC); RowA.clear ();
      }
      RowA.push_back (RowB); RowB.clear ();
    }
    TensorOut.push_back (RowA); RowA.clear ();
  }
  return TensorOut;
}

//------------------------------------------------------------------------------
// updateVoigtMatrix () : Generates the contracted Voigt matrix based on the
// values of the three magnetostriction coefficients. Note that this routine
// will work for crystals of CUBIC SYMMETRY ONLY at the moment. Others must be 
// added to the switch construct if you wish to use them. See Nye.
//
int MgstTensor::updateVoigtMatrix () {
  vector <double> x;  // Matrix element workspace
  int i, j;           // Matrix indicies
  
  // In Voigt matrix notation, the tensor is reduced to a 6x6 matrix.
  // However, the number of independent components for all but monoclinic 
  // structures is even less than 36. 
  
  switch (theSymmetry) {
    
    // For ISOTROPIC symmetry, there are two independent coefficients:
    //
    // / 0 1 1 . . . \ 	All components marked   0 are equal
    // |   0 1 . . . |                          1 are equal
    // |     0 . . . |                          2 equal 2 (#0 - #1)
    // |       2 . . |                          . are zero
    // |         2 . |
    // \           2 /  The matrix is symmetrical about the leading diagonal
    case MT_TENSOR_ISOTROPIC:
      x.push_back ((theLambdaA0 + 2.0 * theLambdaG2) / 3.0);
      x.push_back ((theLambdaA0 - theLambdaG2) / 3.0);
      
      // Fill the top left quadrant of the matrix
      for (i = 0; i < 3; i ++) {
        for (j = i + 1; j < 3; j ++) {
          VoigtMatrix[i][j] = x[1];
          VoigtMatrix[j][i] = x[1];
        }
        VoigtMatrix[i][i] = x[0];
      }
      
      // Fill the rest of the matrix
      for (i = 3; i < 6; i ++) {
        for (j = i + 1; j < 6; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }
        for (j = 0; j < i; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }         
        VoigtMatrix[i][i] = 2.0 * (x[0] - x[1]);
      }
      break;


    // For a CUBIC crystal, the following holds:
    //
    // / 0 1 1 . . . \ 	All components marked   0 are equal
    // |   0 1 . . . |                          1 are equal
    // |     0 . . . |                          2 are equal
    // |       2 . . |                          . are zero
    // |         2 . |
    // \           2 /  The matrix is symmetrical about the leading diagonal
    case MT_TENSOR_CUBIC:
      x.push_back ((theLambdaA0 + 2.0 * theLambdaG2) / 3.0);
      x.push_back ((theLambdaA0 - theLambdaG2) / 3.0);
      x.push_back (2.0 * theLambdaE2);
      
      // Fill the top left quadrant of the matrix
      for (i = 0; i < 3; i ++) {
        for (j = i + 1; j < 3; j ++) {
          VoigtMatrix[i][j] = x[1];
          VoigtMatrix[j][i] = x[1];
        }
        VoigtMatrix[i][i] = x[0];
      }
      
      // Fill the rest of the matrix
      for (i = 3; i < 6; i ++) {
        for (j = i + 1; j < 6; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }
        for (j = 0; j < i; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }         
        VoigtMatrix[i][i] = x[2];
      }
      break;
    
    // The same can be done for CYLINDRICAL symmetry
    //
    // / 0 1 2 . . . |  
    // | 1 . 2 . . . |  Same concept as before, but matrix is not symmetrical
    // | 3 3 4 . . . |
    // | . . . 5 . . | 
    // | . . . . 5 . |
    // | . . . . . 6 /
    
    case MT_TENSOR_CYLINDRICAL:
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        - 1.0 / 9.0 * (theLambda1A2 - 3.0 / 2.0 * theLambda2A2) 
        + 0.5 * theLambdaG2);
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        - 1.0 / 9.0 * (theLambda1A2 - 3.0 / 2.0 * theLambda2A2)
        - 0.5 * theLambdaG2);
      x.push_back (1.0 / 3.0 * (theLambda1A0 + 2.0 * theLambda2A0)
        - 1.0 / 9.0 * (theLambda1A2 - 3.0 * theLambda2A2));
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        + 2.0 / 9.0 * (theLambda1A2 - 3.0 / 2.0 * theLambda2A2));
      x.push_back (1.0 / 3.0 * (theLambda1A0 + 2.0 * theLambda2A0)
        + 2.0 / 9.0 * (theLambda1A2 + 3.0 * theLambda2A2));
      x.push_back (2.0 * theLambdaD2);
      x.push_back (2.0 * (x [0] - x [1]));
      
      // Fill the non-zero elements of the matrix
      VoigtMatrix[0][0] = x[0]; 
      VoigtMatrix[0][1] = x[1]; VoigtMatrix[1][0] = x[1];
      VoigtMatrix[0][2] = x[2]; VoigtMatrix[1][2] = x[2];
      VoigtMatrix[2][0] = x[3]; VoigtMatrix[2][1] = x[3];
      VoigtMatrix[2][2] = x[4];
      VoigtMatrix[3][3] = x[5]; VoigtMatrix[4][4] = x[5];
      VoigtMatrix[5][5] = x[6];
      
      // Fill the rest of the matrix with zeros
      VoigtMatrix[1][1] = 0.0;
      for (i = 3; i < 6; i ++) {
        for (j = i + 1; j < 6; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }
        for (j = 0; j < i; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }         
      }
      break;

    // Now TETRAGONAL symmetry for T1 LAUE GROUPS (422, 4/mmm, 4mm, -42m)
    //
    // / 0 1 2 . . . |  
    // | 1 0 2 . . . |  Similar to cylindrical symmetrical but M22 now used
    // | 3 3 4 . . . |
    // | . . . 5 . . | 
    // | . . . . 5 . |
    // | . . . . . 6 /
    
    case MT_TENSOR_TETRAGONAL_T1:
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        - 1.0 / 6.0 * (theLambda1A2 - theLambda2A2) + 0.5 * theLambdaG2);
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        - 1.0 / 6.0 * (theLambda1A2 - theLambda2A2) - 0.5 * theLambdaG2);
      x.push_back (1.0 / 3.0 * (theLambda1A0 + 2.0 * theLambda2A0)
        - 1.0 / 6.0 * (theLambda1A2 + 2.0 * theLambda2A2));
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        + 1.0 / 3.0 * (theLambda1A2 - theLambda2A2));
      x.push_back (1.0 / 3.0 * (theLambda1A0 + 2.0 * theLambda2A0)
        + 1.0 / 3.0 * (theLambda1A2 + 2.0 * theLambda2A2));
      x.push_back (1.0 / 3.0 * theLambda1A0 + 4.0 / 9.0 * theLambdaE2);
      x.push_back (1.0 / 3.0 * theLambda1A0 + 2.0 / 9.0 * theLambdaD2);
      
      // Fill the non-zero elements of the matrix
      VoigtMatrix[0][0] = x[0]; VoigtMatrix[1][1] = x[0];
      VoigtMatrix[0][1] = x[1]; VoigtMatrix[1][0] = x[1];
      VoigtMatrix[0][2] = x[2]; VoigtMatrix[1][2] = x[2];
      VoigtMatrix[2][0] = x[3]; VoigtMatrix[2][1] = x[3];
      VoigtMatrix[2][2] = x[4];
      VoigtMatrix[3][3] = x[5]; VoigtMatrix[4][4] = x[5];
      VoigtMatrix[5][5] = x[6];
      
      // Fill the rest of the matrix with zeros
      for (i = 3; i < 6; i ++) {
        for (j = i + 1; j < 6; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }
        for (j = 0; j < i; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }         
      }
      break;
    
    case MT_TENSOR_TETRAGONAL_T2:
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        - 1.0 / 6.0 * (theLambda1A2 - theLambda2A2) + 0.5 * theLambdaG2);
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        - 1.0 / 6.0 * (theLambda1A2 - theLambda2A2) - 0.5 * theLambdaG2);
      x.push_back (1.0 / 3.0 * (theLambda1A0 + 2.0 * theLambda2A0)
        - 1.0 / 6.0 * (theLambda1A2 + 2.0 * theLambda2A2));
      x.push_back (1.0 / 3.0 * (theLambda1A0 - theLambda2A0)
        + 1.0 / 3.0 * (theLambda1A2 - theLambda2A2));
      x.push_back (1.0 / 3.0 * (theLambda1A0 + 2.0 * theLambda2A0)
        + 1.0 / 3.0 * (theLambda1A2 + 2.0 * theLambda2A2));
      x.push_back (1.0 / 3.0 * theLambda1A0 + 4.0 / 9.0 * theLambdaE2);
      x.push_back (1.0 / 3.0 * theLambda1A0 + 2.0 / 9.0 * theLambda2G2);
      
      x.push_back (1.0 / 3.0 * theLambda1A0 - 1.0 / 6.0 * theLambda1A2
        + 1.0 / 3.0 * theLambda4G2);
      x.push_back (1.0 / 3.0 * theLambda1A0 - 1.0 / 3.0 * theLambda2A0
        + 1.0 / 3.0 * theLambda3G2);
      
      // Fill the non-zero elements of the matrix
      VoigtMatrix[0][0] = x[0]; VoigtMatrix[1][1] = x[0];
      VoigtMatrix[0][1] = x[1]; VoigtMatrix[1][0] = x[1];
      VoigtMatrix[0][2] = x[2]; VoigtMatrix[1][2] = x[2];
      VoigtMatrix[2][0] = x[3]; VoigtMatrix[2][1] = x[3];
      VoigtMatrix[2][2] = x[4];
      VoigtMatrix[3][3] = x[5]; VoigtMatrix[4][4] = x[5];
      VoigtMatrix[5][5] = x[6];
      
      // Fill the rest of the matrix with zeros
      for (i = 3; i < 6; i ++) {
        for (j = i + 1; j < 6; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }
        for (j = 0; j < i; j ++) {
          VoigtMatrix[i][j] = 0.0;
          VoigtMatrix[j][i] = 0.0;
        }         
      }

      VoigtMatrix[1][6] = x[7]; VoigtMatrix[2][6] = -1.0 * x[7];
      VoigtMatrix[6][1] = x[8]; VoigtMatrix[6][2] = -1.0 * x[8];

      break;
      
    default:
      cout << "Invalid crystal symmetry specified. Aborting!" << endl;
      assert (0);
      return MT_ERROR;
  }
  
  // Add on the values of any Voigt matrix fit paramters
  for (i = 0; i < 6; i ++) {
    for (j = 0; j < 6; j ++) {
      VoigtMatrix [i][j] += VoigtParams [i][j];
    }
  }
  
  // Finally, apply the crystallite preferential orientation to the Voigt matrix
  voigtPrefOrientation ();
  return MT_NO_ERROR;
}

//------------------------------------------------------------------------------
// updateMgstTensor () : Generates the full 4th-rank magnetostriction tensor by
// expanding the contracted Voigt matrix.
//
int MgstTensor::updateMgstTensor (double theVoigtMatrix [6][6]) {
  int i, j;
  double VoigtComponent;

  // Now reconstruct the original tensor from the voigt matrix using the
  // following rules:   (from 'Physical Properties of Crystals', Nye, p134)
  //
  //     s(ijkl) = s(mn)  when m and n are 1, 2, or 3
  //    2s(ijkl) = s(mn)  when either m or n are 4, 5, or 6
  //    4s(ijkl) = s(mn)  when both m and n are 4, 5, or 6
  
  for (i = 0; i < 3; i ++) {
    for (j = 0; j < 3; j ++) {
      Tensor[i][i][j][j] = theVoigtMatrix [i][j];

      VoigtComponent = theVoigtMatrix [(i+2)%3+3][(j+2)%3+3] / 4.0;
      Tensor[i][(i+1)%3][j][(j+1)%3] = VoigtComponent;
      Tensor[(i+1)%3][i][(j+1)%3][j] = VoigtComponent;
      Tensor[i][(i+1)%3][(j+1)%3][j] = VoigtComponent;    
      Tensor[(i+1)%3][i][j][(j+1)%3] = VoigtComponent;

      VoigtComponent = theVoigtMatrix[(i+2)%3+3][j] / 2.0;
      Tensor[i][(i+1)%3][j][j] = VoigtComponent;
      Tensor[(i+1)%3][i][j][j] = VoigtComponent;
      
      VoigtComponent = theVoigtMatrix[i][(j+2)%3+3] / 2.0;
      Tensor[i][i][(j+1)%3][j] = VoigtComponent;
      Tensor[i][i][j][(j+1)%3] = VoigtComponent;
    }
  }
  return MT_NO_ERROR;
}


int MgstTensor::rotate (double x, double y, double z) {
  int i, j, k;
  double Sum;
  double PartialProduct[6];
  double RotatedVoigt[6][6];

  vector < vector <double> > BondMatrix = bondMatrix (eulerMatrix (x, y, z));

  for (i = 0; i < 6; i ++) {
    for (j = 0; j < 6; j ++) {
      Sum = 0.0;
      for (k = 0; k < 6; k ++) {
        Sum += BondMatrix[i][k] * VoigtMatrix [k][j];
      }
      PartialProduct[j] = Sum;
    }
    for (j = 0; j < 6; j ++) {
      Sum = 0.0;
      for (k = 0; k < 6; k ++) {
        Sum += PartialProduct[k] * BondMatrix [j][k];
      }
      RotatedVoigt[i][j] = Sum;
    }
  }
  updateMgstTensor (RotatedVoigt);
  return MT_NO_ERROR;
}

int MgstTensor::voigtPrefOrientation () {
  int i, j, k;
  double Sum;
  double PartialProduct[6];
  double RotatedVoigt[6][6];
    
  vector < vector <double> > BondMatrix = bondMatrix (
    eulerMatrix (PrefOrientation.x, PrefOrientation.y, PrefOrientation.z));

  // Rotate the Voigt matrix using the Bond Matrix
  for (i = 0; i < 6; i ++) {
    for (j = 0; j < 6; j ++) {
      Sum = 0.0;
      for (k = 0; k < 6; k ++) {
        Sum += BondMatrix[i][k] * VoigtMatrix [k][j];
      }
      PartialProduct[j] = Sum;
    }
    for (j = 0; j < 6; j ++) {
      Sum = 0.0;
      for (k = 0; k < 6; k ++) {
        Sum += PartialProduct[k] * BondMatrix [j][k];
      }
      RotatedVoigt[i][j] = Sum;
    }
  }
  
  // This is where things differ from ::rotate in that the original VoigtMatrix
  // is now overwritten so that all future rotations are with respect to the
  // preferential orientation.
  for (i = 0; i < 6; i ++) {
    for (j = 0; j < 6; j ++) {
      VoigtMatrix[i][j] = RotatedVoigt[i][j];
//      cout << VoigtMatrix[i][j] << endl;
    }
  }
  return MT_NO_ERROR;
}  

vector<vector<double> > MgstTensor::eulerMatrix (double x, double y, double z) {
  vector < vector <double> > Matrix;
  vector <double> Row;
  
  Row.push_back (cos(y) * cos(z));
  Row.push_back (cos(y) * sin(z));
  Row.push_back (0.0 - sin(y));
  Matrix.push_back (Row); Row.clear ();
  Row.push_back ((sin(x) * sin(y) * cos(z)) - (cos(x) * sin(z)));
  Row.push_back ((sin(x) * sin(y) * sin(z)) + (cos(x) * cos(z)));
  Row.push_back (sin(x) * cos(y));
  Matrix.push_back (Row); Row.clear ();
  Row.push_back ((cos(x) * sin(y) * cos(z)) + (sin(x) * sin(z)));
  Row.push_back ((cos(x) * sin(y) * sin(z)) - (sin(x) * cos(z)));
  Row.push_back (cos(x) * cos(y));
  Matrix.push_back (Row); Row.clear ();
  
  return Matrix;
}


//------------------------------------------------------------------------------
// bondMatrix (vector < vector <double> >) : Creates the 6x6 Bond transformation
// matrix based on the input 3x3 matrix of direction cosines. These direction
// cosines can be generated with the ::eulerMatrix routine
// 
vector < vector <double> > MgstTensor::bondMatrix (
  vector < vector <double> > DirCosines) {
  
  int i, j;         // Matrix indices
  vector < vector <double> > RtnBondMatrix;
  vector <double> Row;
  double BondMatrix [6][6];
  for (i = 0; i < 3; i ++) {
    for (j = 0; j < 3; j ++) {
      BondMatrix[i][j] = pow (DirCosines[i][j], 2);
      BondMatrix[i][j+3] = DirCosines[i][(j+1)%3] * DirCosines [i][(j+2)%3];
      BondMatrix[i+3][j] = DirCosines[(i+1)%3][j] * DirCosines [(i+2)%3][j]*2.0;
//      a = (i+2)%3; b = (j+2)%3; c = (i+1)%3; d = (i+1)%3;
//      BondMatrix[i+3][j+3] = (DirCosines[a][b] * DirCosines[c][d]) +
//                             (DirCosines[a][d] * DirCosines[c][b]);
    }
  }
  BondMatrix[3][3] 
    = DirCosines[1][1] * DirCosines[2][2] + DirCosines[1][2] * DirCosines[2][1];
  BondMatrix[3][4] 
    = DirCosines[1][0] * DirCosines[2][2] + DirCosines[1][2] * DirCosines[2][0];
  BondMatrix[3][5] 
    = DirCosines[1][1] * DirCosines[2][0] + DirCosines[1][0] * DirCosines[2][1];
  BondMatrix[4][3] 
    = DirCosines[0][1] * DirCosines[2][2] + DirCosines[0][2] * DirCosines[2][1];
  BondMatrix[4][4] 
    = DirCosines[0][0] * DirCosines[2][2] + DirCosines[0][2] * DirCosines[2][0];
  BondMatrix[4][5] 
    = DirCosines[0][0] * DirCosines[2][1] + DirCosines[0][1] * DirCosines[2][0];
  BondMatrix[5][3] 
    = DirCosines[1][1] * DirCosines[0][2] + DirCosines[0][1] * DirCosines[1][2];
  BondMatrix[5][4] 
    = DirCosines[0][0] * DirCosines[1][2] + DirCosines[0][2] * DirCosines[1][0];
  BondMatrix[5][5] 
    = DirCosines[1][1] * DirCosines[0][0] + DirCosines[0][1] * DirCosines[1][0];
  for (i = 0; i < 6; i ++) {
    for (j = 0; j < 6; j ++) {
      Row.push_back (BondMatrix [i][j]);
    }
    RtnBondMatrix.push_back (Row);
    Row.clear ();
  } 
  return RtnBondMatrix;
}


//******************************************************************************
// DEBUG ROUTINES
//******************************************************************************

int MgstTensor::saveVoigt (char *Filename) {
  ofstream VoigtFile (Filename);
  if (! VoigtFile.is_open()) {
    return MT_ERROR;
  }
  
  VoigtFile.precision (MT_OUTPUT_PRECISION);
  for (int m = 0; m < 6; m ++) {
    for (int n = 0; n < 6; n ++) {
      VoigtFile << showpos << scientific << VoigtMatrix [m][n] << " ";
    }
    VoigtFile << endl;
  }
  VoigtFile.close ();
  return MT_NO_ERROR;
}

int MgstTensor::saveTensor (char *Filename) {
  ofstream TensorFile (Filename);
  if (! TensorFile.is_open()) {
    return MT_ERROR;
  }
  
  TensorFile.precision (MT_OUTPUT_PRECISION);
  TensorFile << "i   j   k   l   Tensor Element" << endl;
  for (int i = 0; i < 3; i ++) {
    for (int j = 0; j < 3; j ++) {
      for (int k = 0; k < 3; k ++) {
        for (int l = 0; l < 3; l ++) {
          TensorFile << scientific << i << "   " << j << "   " << k 
            << "   " << l << "   " << Tensor[i][j][k][l] << endl;
        }
      }
    }
  }
  TensorFile.close ();
  return MT_NO_ERROR;
}

int MgstTensor::rotateDebug (double x, double y, double z) {
  saveVoigt ("debug.voigt1");
  saveTensor ("debug.tensor1");

  int m, n;
  vector < vector <double> > EulerMatrix = eulerMatrix (x, y, z);
  vector < vector <double> > BondMatrix = bondMatrix (EulerMatrix);
  rotate (x, y, z);
  
  ofstream FileOut ("debug.euler");
  if (! FileOut.is_open()) {
    return MT_ERROR;
  }
  FileOut.precision (MT_OUTPUT_PRECISION);
  for (m = 0; m < 3; m ++) {
    for (n = 0; n < 3; n ++) {
      FileOut << showpos << scientific << EulerMatrix [m][n] << " ";
    }
    FileOut << endl;
  }
  FileOut.close ();
  
  FileOut.open ("debug.bond");
  if (! FileOut.is_open()) {
    return MT_ERROR;
  }
  FileOut.precision (MT_OUTPUT_PRECISION);
  for (m = 0; m < 6; m ++) {
    for (n = 0; n < 6; n ++) {
      FileOut << showpos << scientific << BondMatrix [m][n] << " ";
    }
    FileOut << endl;
  }
  FileOut.close ();
  
  saveVoigt ("debug.voigt2");
  saveTensor ("debug.tensor2");
  
  return MT_NO_ERROR;
}
