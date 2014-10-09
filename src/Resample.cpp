// Resample.cpp
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
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <vector>
#include <iostream>
#include "DiffXasSpectrum.h"
using namespace::std;

#define SPLINE_ORDER 4  // Cubic spline

//-------------------------------------------------------------------------------
// resample (vector<double> *, vector<double> *, vector<double> *, double) :
// Given a set of coordinates sampled at arbitrary intervals along the x-axis,
// this function uses GSL's B-spline routines to resample the data with points
// spaced at equal intervals. This spacing is defined by arg4, and may be used to
// effectively define the Nyquist frequency for the input signal. All frequencies
// lower than this are representable by the spline and are considered part of the
// signal. All frequencies greater than this are rejected by the spline and so
// may be considered as noise. These residuals are returned in the vector at arg3
//
int resample (vector<double> *XCoords, 
  vector<double> *YCoords, vector<double> *Residuals, double KnotSpacing) {

  vector<double> PaddedXCoords, PaddedYCoords;
  size_t NumDataPoints = XCoords -> size ();
  size_t NumKnots = 
    int ((XCoords -> at (NumDataPoints - 1) - XCoords -> at(0)) / KnotSpacing);
  size_t NumCoefficients = NumKnots + 2;
  /*if (NumDataPoints < NumCoefficients) {
    cout << "ERROR: The experimental spectrum contains too low a density of data points to be resampled prior to noise extraction."
      << endl << "Please increase the sample rate to at least dk=" << KnotSpacing << "(dE=" << KnotSpacing * KnotSpacing * ALPHA << ")."
      << endl << "or recompile DEXA with a lower MAX_RADIUS." << endl;
    return DX_ERROR;
  }*/
  size_t i, j;
  double chisq;
  double MinX = XCoords -> at (0);
  double MaxX = XCoords -> at (NumDataPoints - 1);
  
  // For spectra sampled in E-space, the spacing of data points in k-space will
  // be uneven. At low-K, the sample rate may be too low to be represented by a
  // spline based on MAX_RADIUS. At higher k, the spacing decreases until there
  // is a sufficiently high density of points. Find the transition between these
  // two regions and pad the low density region with extra points based on
  // linear interpolation.
  for (i = 0; i < NumDataPoints - 1; i ++) {
    PaddedXCoords.push_back (XCoords -> at (i));
    PaddedYCoords.push_back (YCoords -> at (i));
    double DeltaX = XCoords -> at (i + 1) - XCoords -> at (i);
    if (DeltaX > KnotSpacing) {
      size_t NumPadPoints = int (DeltaX / KnotSpacing) + 1;
      for (j = 1; j < NumPadPoints; j ++) {
        PaddedXCoords.push_back (XCoords -> at (i) + (DeltaX/NumPadPoints) * j);
        PaddedYCoords.push_back (YCoords -> at (i) + ((double(j) / double(NumPadPoints)) *
          (YCoords -> at (i + 1) - YCoords -> at (i))));
      }
    }
  }

  // Handle the last data point explicitly and reset NumDataPoints so as to 
  // include the newly added padding points  
  PaddedXCoords.push_back (XCoords -> at (NumDataPoints - 1));
  PaddedYCoords.push_back (YCoords -> at (NumDataPoints - 1));
  NumDataPoints = PaddedXCoords.size();

  // Allocate a cubic bspline workspace
  gsl_bspline_workspace *bw = gsl_bspline_alloc(SPLINE_ORDER, NumKnots);
  gsl_vector *B = gsl_vector_alloc(NumCoefficients);
  gsl_vector *x = gsl_vector_alloc(NumDataPoints);
  gsl_vector *y = gsl_vector_alloc(NumDataPoints);
  gsl_vector *c = gsl_vector_alloc(NumCoefficients);
  gsl_matrix *X = gsl_matrix_alloc(NumDataPoints, NumCoefficients);
  gsl_matrix *cov = gsl_matrix_alloc(NumCoefficients, NumCoefficients);
  gsl_multifit_linear_workspace *mw 
    = gsl_multifit_linear_alloc(NumDataPoints, NumCoefficients);
  gsl_vector *Knots = gsl_vector_alloc(NumKnots);

  // Prepare GSL vectors with the input data set and define a knot vector for the
  // resampled data points.
  for (i = 0; i < NumDataPoints; i ++) {
//    cout << PaddedXCoords [i] << ", " << PaddedYCoords[i] << endl;
    gsl_vector_set (x, i, PaddedXCoords[i]);
    gsl_vector_set (y, i, PaddedYCoords[i]);
  }
  for (i = 0; i < NumKnots - 1; i ++) {
    gsl_vector_set (Knots, i, MinX + i * KnotSpacing);
  }
  gsl_vector_set (Knots, NumKnots - 1, MaxX);
  gsl_bspline_knots (Knots, bw);

  // Construct the fit matrix that will be used to optimise the spline
  for (i = 0; i < NumDataPoints; i ++) {
    double xi = gsl_vector_get(x, i);

    /* compute B_j(xi) for all j */
    gsl_bspline_eval(xi, B, bw);

    /* fill in row i of X */
    for (j = 0; j < NumCoefficients; ++j) {
      double Bj = gsl_vector_get(B, j);
      gsl_matrix_set(X, i, j, Bj);
    }
  }

  // Do the fit
  gsl_multifit_linear(X, y, c, cov, &chisq, mw);

  // Get the spectrum residuals i.e. the parts rejected by the fitted spline.
  Residuals -> clear ();
  for (i = 0; i < NumDataPoints; i ++) {
    double YSpline, YError;
    gsl_bspline_eval (PaddedXCoords[i], B, bw);
    gsl_multifit_linear_est (B, c, cov, &YSpline, &YError);
    Residuals -> push_back ((PaddedYCoords[i] - YSpline));
  }

  // Now use the spline to reconstruct the data set at the given sample rate
  XCoords -> clear ();
  YCoords -> clear ();
  double XSpline, YSpline, YError;
  for (XSpline = MinX; XSpline < MaxX; XSpline += KnotSpacing){
    gsl_bspline_eval (XSpline, B, bw);
    gsl_multifit_linear_est (B, c, cov, &YSpline, &YError);
    XCoords -> push_back (XSpline);
    YCoords -> push_back (YSpline);
  }

  // Handle the final data point explicitly
  XSpline = MaxX;
  gsl_bspline_eval (XSpline, B, bw);
  gsl_multifit_linear_est (B, c, cov, &YSpline, &YError);
  XCoords -> push_back (XSpline);
  YCoords -> push_back (YSpline);

  // Clean up the GSL workspace and exit
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);
  gsl_vector_free(Knots);
  return DX_NO_ERROR;
}


//-------------------------------------------------------------------------------
// processNoise(...) : This function performs a similar operation to resample()
// above, but works on the noise spectrum ie the residuals from resample (). The
// noise's key difference is that a bandwidth correction is first applied to the
// noise based upon any Fourier filtering performed on the original signal from
// which the noise came. This accounts for the noise at frequencies lower than
// the Nyquist limit for the signal, but greater than the upper limit of the
// Fourier window.
//
void processNoise (vector<double> *XCoords, vector<double> *YCoords, 
  double FilterMaxR, double AbsMaxR, vector<double> *SamplePoints) {
  
  vector<double> PaddedXCoords, PaddedYCoords;
  double KnotSpacing = SamplePoints -> at(1) - SamplePoints -> at(0);
  unsigned int NumDataPoints = XCoords -> size ();
  unsigned int NumKnots = SamplePoints -> size ();
  unsigned int NumCoefficients = NumKnots + 2;
  double chisq;
    
  // Prepare the residuals. First apply a bandwidth correction to account for
  // the noise in the part of the spectrum between the upper-R limit of the
  // Fourier filter and the low-R limit of the noise extraction in resample().
//  for (unsigned int i = 0; i < NumDataPoints; i ++) {
//    YCoords -> at (i) = /*pow(*/YCoords -> at(i) * FilterMaxR / AbsMaxR/*, 2.0)*/;
//  }

  // For spectra sampled in E-space, the spacing of data points in k-space will
  // be uneven. At low-K, the sample rate may be too low to be represented by a
  // spline based on MAX_RADIUS. At higher k, the spacing decreases until there
  // is a sufficiently high density of points. Find the transition between these
  // two regions and pad the low density region with extra points based on
  // linear interpolation.
  for (unsigned int i = 0; i < NumDataPoints - 1; i ++) {
    PaddedXCoords.push_back (XCoords -> at (i));
    PaddedYCoords.push_back (YCoords -> at (i));
    double DeltaX = XCoords -> at (i + 1) - XCoords -> at (i);
    if (DeltaX > KnotSpacing) {
      size_t NumPadPoints = int (DeltaX / KnotSpacing) + 1;
      for (unsigned int j = 1; j < NumPadPoints; j ++) {
        PaddedXCoords.push_back (XCoords -> at (i) + (DeltaX/NumPadPoints) * j);
        PaddedYCoords.push_back (YCoords -> at (i) + ((double(j) / double(NumPadPoints)) *
          (YCoords -> at (i + 1) - YCoords -> at (i))));
      }
    }
  }

  // Handle the last data point explicitly and reset NumDataPoints so as to 
  // include the newly added padding points  
  PaddedXCoords.push_back (XCoords -> at (NumDataPoints - 1));
  PaddedYCoords.push_back (YCoords -> at (NumDataPoints - 1));
  NumDataPoints = PaddedXCoords.size();


  // allocate a cubic bspline workspace
  gsl_bspline_workspace *bw = gsl_bspline_alloc(SPLINE_ORDER, NumKnots);
  gsl_vector *B = gsl_vector_alloc(NumCoefficients);
  gsl_vector *x = gsl_vector_alloc(NumDataPoints);
  gsl_vector *y = gsl_vector_alloc(NumDataPoints);
  gsl_vector *c = gsl_vector_alloc(NumCoefficients);
  gsl_matrix *X = gsl_matrix_alloc(NumDataPoints, NumCoefficients);
  gsl_matrix *cov = gsl_matrix_alloc(NumCoefficients, NumCoefficients);
  gsl_multifit_linear_workspace *mw
    = gsl_multifit_linear_alloc(NumDataPoints, NumCoefficients);
  gsl_vector *Knots = gsl_vector_alloc(NumKnots);

  /* this is the data to be fitted */
  for (unsigned int i = 0; i < NumDataPoints; i ++) {
    gsl_vector_set (x, i, PaddedXCoords[i]);
    gsl_vector_set (y, i, PaddedYCoords[i]);
  }
  for (unsigned int i = 0; i < NumKnots; i ++) {
    gsl_vector_set (Knots, i, SamplePoints -> at (i));
  }
  gsl_bspline_knots (Knots, bw);
  /* construct the fit matrix X */
  for (unsigned int i = 0; i < NumDataPoints; i ++) {
    double xi = gsl_vector_get(x, i);

    /* compute B_j(xi) for all j */
    gsl_bspline_eval(xi, B, bw);

    /* fill in row i of X */
    for (unsigned int j = 0; j < NumCoefficients; ++j) {
      double Bj = gsl_vector_get(B, j);
      gsl_matrix_set(X, i, j, Bj);
    }
  }
  // Do the fit
  gsl_multifit_linear(X, y, c, cov, &chisq, mw);

  // Now use the spline to reconstruct the data set at the given sample rate
  XCoords -> clear ();
  YCoords -> clear ();
  for (unsigned int i = 0; i < NumKnots; i ++){
    double XSpline, YSpline, YError;
    XSpline = SamplePoints -> at (i);
    gsl_bspline_eval (XSpline, B, bw);
    gsl_multifit_linear_est (B, c, cov, &YSpline, &YError);
    XCoords -> push_back (XSpline);
    YCoords -> push_back (YSpline * GSL_SIGN(YSpline));
  }

  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);
  gsl_vector_free(Knots);

}








