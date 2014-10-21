// ThermalSpectrum.h
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
#ifndef THERMALSPECTRUM_H_
#define THERMALSPECTRUM_H_

#include "PathFinder.h"
#include "DiffXasSpectrum.h"

#define TH_NO_ERROR 0
#define TH_ERROR    1


// Important: Beamline frame of reference is
//    +z beam propagation direction,
//    +y from hutch floor vertically up towards ceiling
//    +x away from the storage ring in direction orthogonal to +z and +y

class ThermalSpectrum :
  public DiffXasSpectrum,
  public PathFinder {

public:
	ThermalSpectrum () { SpectrumReady = false; dChiCheb = gsl_cheb_alloc (CHEB_ORDER); }
  ~ThermalSpectrum () { /*gsl_cheb_free (dChiCheb);*/ }
  int copySpectrum (ThermalSpectrum Input);

  // Local declarations of DiffXasSpectrum abstract base functions
  double dChi (double x);
  FFData getDChiPath (int PathNo);

private:
  void calculateDChi ();
  vector <Coordinate> doDChiCalculation (vector <Path> paths, vector <Coordinate> XCoords);
  Coordinate doDChiCalculation1 (vector <Path> paths, double k);
  static double dChiChebFunction (double k, void *p);
  void generateDChiCheb (vector <Coordinate> &Data);

  gsl_cheb_series *dChiCheb;  // A Chebyshev approximation of dChi(k)
};

#endif /* THERMALSPECTRUM_H_ */
