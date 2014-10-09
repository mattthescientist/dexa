#include <iostream>
/* deriv/deriv.c
 *
 * Copyright (C) 2004 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

typedef struct parameter_function_struct {
  vector <double> (* function) (double x, vector<double> y, void * params);
  void * params;
} parameter_function;

#define PARAM_FN_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)


static void parameter_deriv (const parameter_function * f, vector<double> x, double h,
  vector<double> *result, vector<double> *abserr_round, 
  vector<double> *abserr_trunc) {

  // Compute the derivative using the 5-point rule (x-h, x-h/2, x,
  // x+h/2, x+h). Note that the central point is not used.
  // Compute the error using the difference between the 5-point and
  // the 3-point rule (x-h,x,x+h). Again the central point is not
  // used.

  vector<double> fm1 = PARAM_FN_EVAL (f, 0.0 - h, x);
  vector<double> fp1 = PARAM_FN_EVAL (f, 0.0 + h, x);

  vector<double> fmh = PARAM_FN_EVAL (f, 0.0 - h / 2, x);
  vector<double> fph = PARAM_FN_EVAL (f, 0.0 + h / 2, x);

  vector<double> r3, r5, e3, e5;
  for (unsigned int i = 0; i < x.size (); i ++) {
    r3.push_back (0.5 * (fp1[i] - fm1[i]));
    r5.push_back ((4.0 / 3.0) * (fph[i] - fmh[i]) - (1.0 / 3.0) * r3[i]);

    e3.push_back ((fabs (fp1[i]) + fabs (fm1[i])) * GSL_DBL_EPSILON);
    e5.push_back (2.0 * (fabs (fph[i]) + fabs (fmh[i])) * GSL_DBL_EPSILON+e3[i]);

    // The truncation error in the r5 approximation itself is O(h^4).
    // However, for safety, we estimate the error from r5-r3, which is
    // O(h^2). By scaling h we will minimise this estimated error, not
    // the actual truncation error in r5.

    result -> push_back(r5[i] / h);
    abserr_trunc -> push_back(fabs ((r5[i] - r3[i]) / h));
    abserr_round -> push_back(fabs (e5[i] / h));
  }
}


int deriv_parameter_central (const parameter_function * f, vector<double> x, 
  double h, vector<double> *result, vector<double> *abserr) {
  
  vector<double> r_0, round, trunc, error, x_processing;
//  vector<double *> r_0_processing, error_processing;
  parameter_deriv (f, x, h, &r_0, &round, &trunc);
  for (unsigned int i = 0; i < x.size (); i ++) {
    error.push_back(round[i] + trunc[i]);
/*    if (round[i] < trunc[i] && (round[i] > 0 && trunc[i] > 0)) {
      x_processing.push_back (x[i]);
      r_0_processing.push_back (&r_0[i]);
      error_processing.push_back (&error[i]);
    }*/
  }

  if (round[0] < trunc[0] && (round[0] > 0 && trunc[0] > 0)) {
    vector<double> r_opt, round_opt, trunc_opt, error_opt;

    // Compute an optimised stepsize to minimize the total error,
    // using the scaling of the truncation error (O(h^2)) and
    // rounding error (O(1/h)).

    double h_opt = h * pow (round[0] / (2.0 * trunc[0]), 1.0 / 3.0);
    parameter_deriv (f, x/*_processing*/, h_opt, &r_opt, &round_opt, &trunc_opt);
    for (unsigned int i = 0; i < x/*_processing*/.size (); i ++) {
      error_opt.push_back(round_opt[i] + trunc_opt[i]);

      // Check that the new error is smaller, and that the new derivative
      // is consistent with the error bounds of the original estimate.

      if (error_opt[i] < /***/error/*_processing*/[i] 
        && fabs (r_opt[i] - /***/r_0/*_processing*/[i]) < 4.0 * /***/error/*_processing*/[i]) {
        /***/r_0/*_processing*/[i] = r_opt[i];
        /***/error/*_processing*/[i] = error_opt[i];
      }
    }
  }
  for (unsigned int i = 0; i < x.size (); i ++) {
    result -> push_back(r_0[i]);
    abserr -> push_back(error[i]);
  }

  return GSL_SUCCESS;
} 



