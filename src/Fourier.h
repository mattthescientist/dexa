// Fourier.h
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
#ifndef FOURIERH
#define FOURIERH

#include "DiffXasSpectrum.h"

#ifndef PI
  #define PI 3.1415926535897932384626433832795
#endif

#ifndef DISABLED
  #define DISABLED -1
#endif

// In which direction is the FFT being performed?
#define FORWARD 1
#define BACKWARD -1

// NUM_FT_POINTS defines the number of data points in the FFT. This MUST be an
// integer power of 2, the value of which is stored in FT_LOG_BASE. In other
// words the following expression MUST hold:  
//
//   2 ^ FT_LOG_BASE = NUM_FT_POINTS;
//
#define NUM_FT_POINTS 4096
#define FT_LOG_BASE 12

class Fourier {
  public:
  
    Fourier ();   // Default constructor
    ~Fourier ();  // Default destructor
    
    void transform (vector<Coordinate> &DataIn, int Direction, double MinK, vector<Coordinate> &DataOut);
    void getTotalFT (vector<Coordinate> &ComplexSpectrum, vector<Coordinate> &TotalFT);
    int hannWindow (vector <double> kMin, vector <double> kMax, vector <double> dK, vector <Coordinate> &Data, const char *OutputName = "");
    int hannWindowComplex (vector <double> kMin, vector <double> kMax, vector <double> dK, vector <Coordinate> &Data, const char *OutputName = "");

  private:
    int saveData (const char *Filename, vector<Coordinate> Data);
};

#endif // FOURIERH
