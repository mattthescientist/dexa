// Fourier.cpp
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
//******************************************************************************
//
#include "FFTRealFixLen.h"
#include "Fourier.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace::std;

////////////////////////////////////////////////////////////////////////////////
// Default Fourier constructor
//
Fourier::Fourier () {
  // Nothing to do
}


////////////////////////////////////////////////////////////////////////////////
// Default Fourier destructor
//
Fourier::~Fourier () {
  // Nothing to do
}


////////////////////////////////////////////////////////////////////////////////
// transform
// Takes the input vector<Coordinate> at arg1 and calculates either its forward
// or backward Fourier transform, as determined by Direction at arg2. The result
// is returned to the calling function in DataOut at arg4. MinK at arg3 is
// needed for the backward transform only, and should be given the value of the
// x ordinate of the first data point in the k-space spectrum (i.e. the one
// initially used to generate the R-space data). This ensures the q-space data 
// is offset correctly from zero when the k-space spectrum starts at k not equal
// to zero. The FFT itself is performed by the FFTReal v2 library by Laurent de 
// Soras (http://desoras.free.fr/prod.html)
//
void Fourier::transform (vector<Coordinate> &DataIn, int Direction, double MinK, vector<Coordinate> &DataOut) {
  unsigned int i;
  Coordinate CurrentPoint;
  float RawData [NUM_FT_POINTS];
  float FTData [NUM_FT_POINTS];
  double DeltaK = DataIn [1].x - DataIn [0].x;

  // First prepare the input FFT Data array by extracting the absorption data
  // from the input spectrum, and padding it with extra zeros (if required) in
  // order to ensure the number of data points is an integer power of 2.
  
  for (i = 0; i < DataIn.size (); i ++) {
    RawData [i] = float(DataIn [i].y);
  }
  for (i = DataIn.size (); i < NUM_FT_POINTS; i ++) {
    RawData [i] = 0.0;
  }
  
  // Now do the Fourier Transform using the FFTReal library
  FFTRealFixLen <FT_LOG_BASE> FFT;
  DataOut.clear ();
  if (Direction == FORWARD) {
    // Forward FFT
    FFT.do_fft (FTData, RawData);
    for (i = 0; i < NUM_FT_POINTS; i ++) {
      CurrentPoint.x = i * PI * (1.0 / (NUM_FT_POINTS * DeltaK));
      CurrentPoint.y = double(FTData [i]);
      DataOut.push_back(CurrentPoint);
    }
  } else {
    // Inverse FFT
    FFT.do_ifft (RawData, FTData);
    FFT.rescale (FTData);
    for (i = 0; i < NUM_FT_POINTS; i ++) {
      CurrentPoint.x = i * PI * (1.0 / (NUM_FT_POINTS * DeltaK)) + MinK;
      CurrentPoint.y = double(FTData [i]);
      DataOut.push_back(CurrentPoint);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// getTotalFT
// Extracts the total Fourier transform from the real and imaginary transform
// components stored in an FFTReal complex array. This is needed to obtain the
// R-space data.
void Fourier::getTotalFT (vector<Coordinate> &ComplexFT,
  vector<Coordinate> &TotalFT) {
  unsigned int i;
  unsigned int NumFTPoints = ComplexFT.size();
  Coordinate CurrentPoint;

  TotalFT.clear();
  for (i = 0; i < NumFTPoints / 2; i ++) {
    CurrentPoint = ComplexFT [i];
    CurrentPoint.y = sqrt (
      pow (CurrentPoint.y, 2) + pow (ComplexFT[(NumFTPoints / 2) + i].y, 2));
    TotalFT.push_back (CurrentPoint);
  }
}


//------------------------------------------------------------------------------
// hannWindow (double, double, double, vector <Coordinate> &) : This function
// applies a Hann window to the input (x, y) coordinate data set passed in at
// arg4. This function works on real-space data for the forward FT from k- to
// R-space. The y ordinate becomes:
//
//  y *= 0.0                                  For                x < kMin - dK/2
//  y *= 0.5 (1 + sin((x - kMin)/ dK * Pi))   For kMin - dK/2 <= x < kMin + dK/2
//  y *= 1.0                                  For kMin + dK/2 <= x < kMax - dK/2
//  y *= 0.5 (1 - sin((x - kMin)/ dK * Pi))   For kMax - dK/2 <= x < kMax + dK/2
//  y *= 0.0                                  For kMax + dK/2 <= x
//
int Fourier::hannWindow 
  (vector <double> kMin, vector <double> kMax, vector <double> dK, vector <Coordinate> &Data, const char *OutputName) {
  unsigned int i, j;
  
  // First initalise the window. The amplitude at all points is set to zero.
  vector<Coordinate> Window;
  Coordinate NewPoint;
  NewPoint.y = 0.0;
  for (i = 0; i < Data.size (); i ++) {
    NewPoint.x = Data[i].x;
    Window.push_back (NewPoint);
  }
  
  // Now cycle through each of the windows specified in the input vectors and
  // adjust the window amplitude accordingly. Don't worry about the amplitude
  // of overlapping windows summing to more than 1 just yet.
  for (j = 0; j < kMin.size (); j ++) {
    if (kMax[j] - dK[j] * 0.5 < kMin[j] + dK[j] * 0.5) {
      cout << "Error: kMax - dK/2 is less than kMin + dk/2" << endl;
      return DX_ERROR;
    }
    if (dK[j] < 0.0) { 
      cout << "Error: dK must be greater than or equal to zero" << endl;
      return DX_ERROR;
    }
    
    for (i = 0; Data[i].x < kMin[j] - (dK[j] * 0.5) && i < Data.size(); i ++) {}
    for (i = i; Data[i].x < kMin[j] + (dK[j] * 0.5) && i < Data.size(); i ++) {
      Window[i].y += 0.5 * (1.0 + sin ((Data[i].x - kMin[j]) / dK[j] * PI));
      }
    for (i = i; Data[i].x < kMax[j] - (dK[j] * 0.5) && i < Data.size(); i ++) {
      Window[i].y += 1.0;
    }
    for (i = i; Data[i].x < kMax[j] + (dK[j] * 0.5) && i < Data.size(); i ++) {
      Window[i].y += 0.5 * (1.0 - sin ((Data[i].x - kMax[j]) / dK[j] * PI));
    }
  }
  
  // Finally, check that the window amplitude never exceeds 1 (as it could if
  // two windows overlap), and apply the window to the data by multiplying its
  // amplitude by that of the window.
  for (i = 0; i < Data.size (); i ++) {
    if (Window[i].y > 1.0) { Window[i].y = 1.0; }
    Data[i].y *= Window[i].y;
  }
  if (OutputName != "") saveData (OutputName, Window);
  return DX_NO_ERROR;
}

//------------------------------------------------------------------------------
// hannWindowComplex (double, double, double, vector <Coordinate> &) : This 
// function applies a Hann window to the input (x, y) coordinate data set passed 
// in at arg4. This function works on the FFTReal complex data in R-space. It is
// therefore only used for the backward FT from R- to q-space.
// The y ordinate becomes:
//
//  y *= 0.0                                  For                x < kMin - dK/2
//  y *= 0.5 (1 + sin((x - kMin)/ dK * Pi))   For kMin - dK/2 <= x < kMin + dK/2
//  y *= 1.0                                  For kMin + dK/2 <= x < kMax - dK/2
//  y *= 0.5 (1 - sin((x - kMin)/ dK * Pi))   For kMax - dK/2 <= x < kMax + dK/2
//  y *= 0.0                                  For kMax + dK/2 <= x
//
int Fourier::hannWindowComplex
  (vector <double> kMin, vector <double> kMax, vector <double> dK, vector <Coordinate> &Data, const char *OutputName) {
  unsigned int i, j;
  unsigned int Len = Data.size () / 2;
  
  // First initalise the window. The amplitude at all points is set to zero.
  vector<Coordinate> Window;
  Coordinate NewPoint;
  NewPoint.y = 0.0;
  for (i = 0; i < Len; i ++) {
    NewPoint.x = Data[i].x;
    Window.push_back (NewPoint);
  }
  
  // Now cycle through each of the windows specified in the input vectors and
  // adjust the window amplitude accordingly. Don't worry about the amplitude
  // of overlapping windows summing to more than 1 just yet.
  for (j = 0; j < kMin.size (); j ++) {
    if (kMax[j] - dK[j] * 0.5 < kMin[j] + dK[j] * 0.5) {
      cout << "Error: kMax - dK/2 is less than kMin + dk/2" << endl;
      return DX_ERROR;
    }
    if (dK[j] < 0.0) { 
      cout << "Error: dK must be greater than or equal to zero" << endl;
      return DX_ERROR;
    }

    for (i = 0; Data[i].x < kMin[j] - (dK[j] * 0.5) && i < Len; i ++) {}
    for (i = i; Data[i].x < kMin[j] + (dK[j] * 0.5) && i < Len; i ++) {
      Window[i].y += 0.5 * (1.0 + sin ((Data[i].x - kMin[j]) / dK[j] * PI));
      }
    for (i = i; Data[i].x < kMax[j] - (dK[j] * 0.5) && i < Len; i ++) {
      Window[i].y += 1.0;
    }
    for (i = i; Data[i].x < kMax[j] + (dK[j] * 0.5) && i < Len; i ++) {
      Window[i].y += 0.5 * (1.0 - sin ((Data[i].x - kMax[j]) / dK[j] * PI));
    }
  }
  
  // Finally, check that the window amplitude never exceeds 1 (as it could if
  // two windows overlap), and apply the window to the data by multiplying its
  // amplitude by that of the window.
  for (i = 0; i < Len; i ++) {
    if (Window[i].y > 1.0) { Window[i].y = 1.0; }
    Data[i].y *= Window[i].y;
    Data[Len + i].y *= Window[i].y;
  }
  if (OutputName != "") saveData (OutputName, Window);
  return DX_NO_ERROR;
}





int Fourier::saveData (const char *Filename, vector<Coordinate> Data) {

  ofstream SpectrumFile (Filename, ios::out);
  if (! SpectrumFile.is_open()) {
    return DX_ERROR;
  }
  
  SpectrumFile.precision (5);
  for (unsigned int i = 0; i < Data.size (); i ++) {
    SpectrumFile << showpos << scientific 
      << Data[i].x << '\t' 
      << Data[i].y << '\t' 
      << endl;
  }
  SpectrumFile.close();
  return DX_NO_ERROR;
}
