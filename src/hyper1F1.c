/*! \file
  \brief Confluent hypergeometric function 1F1(a, b, z) for a = 1 and complex b & z.
  Both Kummer series and continued fractions are used.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sum.h>

#include "hyper1F1.h"

#ifndef CMPLX
#define CMPLX(x, y) ((double complex)((x) + _Complex_I * (y)))
#endif

struct quad_params {
  complex_double b;
  complex_double z;
  int part;
};

/*******************************************************************/

double min(double x1, double x2)
{
  return (x1 < x2) ? x1 : x2;
}

/*******************************************************************/

double exp1mt(double t, void *params)
{
  struct quad_params *qp = (struct quad_params *) params;

  complex_double ans = (qp->b - 1.0) * cexp(qp->z * t) * cpow(1.0 - t, qp->b - 2.0);

  if (qp->part == 0) {
    return creal(ans);
  }
  return cimag(ans);
}

/*******************************************************************/

void hypergeometric1f1_quad(double *b_re, double *b_im,
                            double *z_re, double *z_im,
                            double *f_re, double *f_im)
{
  // computes function 1F1(a,b,z) for a = 1 and complex b & z by quadrature
  // must be optimized: avoid memory allocation!

  gsl_set_error_handler_off();

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im);

  size_t limit = 100;
  double epsabs = 1.0e-12, epsrel = 1.0e-12, err;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);

  struct quad_params qp = {b, z, 0};

  gsl_function F;
  F.function = &exp1mt;
  F.params = &qp;

  qp.part = 0;
  gsl_integration_qag(&F, 0.0, 1.0, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, f_re, &err);

  qp.part = 1;
  gsl_integration_qag(&F, 0.0, 1.0, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, f_im, &err);

  gsl_integration_workspace_free(w);
}

/*******************************************************************/

void hypergeometric1f1_kummer_nmax(double *b_re, double *b_im,
                                   double *z_re, double *z_im,
                                   double *f_re, double *f_im)
{
  // computes function 1F1(a,b,z) for a = 1 and complex b & z by Kummer series

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im);

  int N = (int) ceil(-20.0 / log10(cabs(z / b))) + 5;

  // fprintf (stdout, "\nhypergeometric1f1_kummer: N = %d", N);

  if (N < 1 || N > 1.0e6) {
    fprintf(stdout, "\nwarning: hypergeometric1F1_kummer: N=%d", N);
    fprintf(stdout, "\nb=%le %le\nz=%le %le\nabs(z/b)=%le",
            creal(b), cimag(b), creal(z), cimag(z), cabs(z / b));
  }

  complex_double term = z / (b + N);

  for (int n = N - 1; n >= 0; n--) {
    term = 1.0 + z / (b + n) * term;
  }

  *f_re = creal(term);
  *f_im = cimag(term);
}

/*******************************************************************/

void hypergeometric1f1_kummer_ada(double *b_re, double *b_im,
                                  double *z_re, double *z_im,
                                  double *f_re, double *f_im)
{
  // computes function 1F1(a,b,z) for a = 1 and complex b & z by Kummer series

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im);

  complex_double term, S1, S2;
  double err, eps = DBL_EPSILON;
  int n, Nmax, maxNmax = 1.0e8;

  for (Nmax = 4; Nmax < maxNmax; Nmax *= 2) {
    n = Nmax;
    term = z / (b + n);

    for (n = Nmax - 1; n >= 0; n--) {
      term = 1.0 + z / (b + n) * term;
    }

    S1 = term;

    n = Nmax + 1;
    term = z / (b + n);

    for (n = Nmax; n >= 0; n--) {
      term = 1.0 + z / (b + n) * term;
    }

    S2 = term;

    err = min(cabs((S2 - S1) / S2), cabs(S2 - S1));
    if (err < eps) {
      break;
    }
  }

  // fprintf (stdout, "\nhypergeometric1f1_kummer_ada: "
  // "Nmax = %d err = %le\nS2 = %.20le %.20le", Nmax, err, creal(S2), cimag(S2));

  if (Nmax >= maxNmax) {
    fprintf(stderr,
            "\nwarning: hypergeometric1f1_kummer_ada: "
            "Nmax = %d err = %le\nS2 = %le %le",
            Nmax, err, creal(S2), cimag(S2));
    fprintf(stderr, "\nb = %le %le z = %le %le",
            creal(b), cimag(b), creal(z), cimag(z));
  }

  *f_re = creal(S2);
  *f_im = cimag(S2);
}

/*******************************************************************/

void hypergeometric1f1_kummer_modified_1(double *b_re, double *b_im,
                                         double *z_re, double *z_im,
                                         double *f_re, double *f_im)
{
  // computes modified function 1F1m(a,b,z) for a = 1 and complex b & z
  // by Kummer series 1F1 = 1 + z/b + z^2/b/(b+1)*1F1m

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im);

  int N = (int) ceil(-20.0 / log10(cabs(z / b))) + 5;

  if (N < 1 || N > 1.0e6) {
    fprintf(stdout, "\nwarning: hypergeometric1F1_kummer: N=%d", N);
    fprintf(stdout, "\nb=%le %le\nz=%le %le\nabs(z/b)=%le",
            creal(b), cimag(b), creal(z), cimag(z), cabs(z / b));
  }

  complex_double term = z / (b + N);

  for (int n = N - 1; n >= 2; n--) {
    term = 1.0 + z / (b + n) * term;
  }

  *f_re = creal(term);
  *f_im = cimag(term);
}

/*******************************************************************/

void hypergeometric1f1_kummer_modified_0_nmax(double *b_re, double *b_im,
                                              double *z_re, double *z_im,
                                              double *f_re, double *f_im)
{
  // computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by
  // kummer series 1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im);

  int N = (int) ceil(-20.0 / log10(cabs(z / b))) + 5;

  if (N < 1 || N > 1.0e6) {
    fprintf(stdout, "\nwarning: hypergeometric1F1_kummer: N=%d", N);
    fprintf(stdout, "\nb=%le %le\nz=%le %le\nabs(z/b)=%le",
            creal(b), cimag(b), creal(z), cimag(z), cabs(z / b));
  }

  complex_double term = z / (b + N);

  for (int n = N - 1; n > 2; n--) {
    term = 1.0 + z / (b + n) * term;
  }

  term *= z / (b + 2.0);

  *f_re = creal(term);
  *f_im = cimag(term);
}

/*******************************************************************/

void hypergeometric1f1_kummer_modified_0_ada(double *b_re, double *b_im,
                                             double *z_re, double *z_im,
                                             double *f_re, double *f_im)
{
  // computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by
  // kummer series 1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im);

  complex_double term, S1, S2;
  double err, eps = DBL_EPSILON;
  int n, Nmax, maxNmax = 1.0e8;

  for (Nmax = 4; Nmax < maxNmax; Nmax *= 2) {
    n = Nmax;
    term = z / (b + n);

    for (n = Nmax - 1; n > 2; n--) {
      term = 1.0 + (z / (b + n)) * term;
    }

    S1 = term * (z / (b + 2.0));

    n = Nmax + 1;
    term = z / (b + n);

    for (n = Nmax; n > 2; n--) {
      term = 1.0 + (z / (b + n)) * term;
    }

    S2 = term * (z / (b + 2.0));

    // err = min(abs((S2 - S1) / S2), abs(S2 - S1));
    err = cabs((S2 - S1) / S2);

    if (err < eps) {
      break;
    }
  }

  if (Nmax >= maxNmax) {
    fprintf(stderr,
            "\nwarning: hypergeometric1f1_kummer_modified_0_ada: Nmax = %d "
            "err = %le\nS2 = %le %le",
            Nmax, err, creal(S2), cimag(S2));
    fprintf(stderr, "\nb = %le %le z = %le %le",
            creal(b), cimag(b), creal(z), cimag(z));
  }

  *f_re = creal(S2);
  *f_im = cimag(S2);
}

/*******************************************************************/

void hypergeometric1f1_kummer_modified_0_accel(double *b_re, double *b_im,
                                               double *z_re, double *z_im,
                                               double *f_re, double *f_im)
{
  // computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by
  // Kummer series 1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)
  // do not use: gives unstable results with wrong error estimation!!!

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im);

  int N = (int) min(50, ceil(-20.0 / log10(cabs(z / b))) + 5);

  complex_double term[N];
  double t_re[N], t_im[N];

  int n = 0;
  term[n] = z / (b + 3.0);
  t_re[n] = creal(term[n]);
  t_im[n] = cimag(term[n]);

  for (n = 1; n < N; n++) {
    term[n] = term[n - 1] * (z / (b + n + 3));
    t_re[n] = creal(term[n]);
    t_im[n] = cimag(term[n]);
  }

  gsl_sum_levin_u_workspace *w = gsl_sum_levin_u_alloc((size_t) N);

  double err;

  gsl_sum_levin_u_accel(t_re, (size_t) N, w, f_re, &err);
  if (err > 1.0e-16) {
    fprintf(stdout, "\nerr = %.16le sum_re = %.16le using %ld terms",
            err, *f_re, w->terms_used);
  }

  gsl_sum_levin_u_accel(t_im, (size_t) N, w, f_im, &err);
  if (err > 1.0e-16) {
    fprintf(stdout, "\nerr = %.16le sum_im = %.16le using %ld terms",
            err, *f_im, w->terms_used);
  }

  gsl_sum_levin_u_free(w);
}

/*******************************************************************/

void hypergeometric1f1_cont_fract_1_modified_0_ada(double *b_re, double *b_im,
                                                   double *z_re, double *z_im,
                                                   double *f_re, double *f_im)
{
  // computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by
  // continued fraction 1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)

  complex_double b = CMPLX(*b_re, *b_im), z = CMPLX(*z_re, *z_im), F11m = CMPLX(0.0, 0.0);

  if (cabs(z / b) < 0.1e0) {
    hypergeometric1f1_kummer_modified_0_ada(b_re, b_im, z_re, z_im, f_re, f_im);
  } else {
    // big numbers substraction - better to implement direct continued fraction!
    hypergeometric1f1_cont_fract_1_inv_ada(b_re, b_im, z_re, z_im, f_re, f_im);
    F11m = (*f_re + *f_im * I - 1.0 - z / b) * (b / z) * ((b + 1.0) / z) - 1.0;
    *f_re = creal(F11m);
    *f_im = cimag(F11m);
  }
}

/*******************************************************************/

void hypergeometric1f1_cont_fract_1_inv_ada(double *b_re, double *b_im,
                                            double *z_re, double *z_im,
                                            double *f_re, double *f_im)
{
  // computes function 1F1(a,b,z) for a = 1 and complex b & z by inverse
  // evaluation of continued fractions
  complex_double b = CMPLX(*b_re - 1.0, *b_im), z = CMPLX(*z_re, *z_im);

  complex_double term, S1, S2;
  double err, eps = DBL_EPSILON;
  int n, Nmax, maxNmax = 1.0e6;

  for (Nmax = 4; Nmax < maxNmax; Nmax *= 2) {
    n = Nmax;
    term = n * z / (b - z + n);
    for (n = Nmax - 1; n > 0; n--) {
      term = n * z / (b - z + n + term);
    }

    S1 = b / (b - z + term);

    n = Nmax + 1;
    term = n * z / (b - z + n);
    for (n = Nmax; n > 0; n--) {
      term = n * z / (b - z + n + term);
    }

    S2 = b / (b - z + term);

    // err = min(cabs((S2 - S1) / S2), cabs(S2 - S1));
    err = cabs((S2 - S1) / S2);

    if (err < eps) {
      break;
    }
  }

  // fprintf(stdout, "\nhypergeometric1f1_cont_fract_1_inv_ada: "
  // "Nmax = %d err = %le\nS2 = %.20le %.20le", Nmax, err, creal(S2), cimag(S2));

  if (Nmax >= maxNmax) {
    fprintf(stderr,
            "\nwarning: hypergeometric1f1_cont_fract_1_inv_ada: Nmax = %d err "
            "= %le\nS2 = %le %le",
            Nmax, err, creal(S2), cimag(S2));
    fprintf(stderr, "\nb = %le %le z = %le %le term = %le %le",
            creal(b) + 1.0, cimag(b), creal(z), cimag(z), creal(term), cimag(term));
  }

  *f_re = creal(S2);
  *f_im = cimag(S2);
}

/*******************************************************************/

void hypergeometric1f1_cont_fract_1_inv_nmax(double *b_re, double *b_im,
                                             double *z_re, double *z_im,
                                             double *f_re, double *f_im)
{
  // computes function 1F1(a,b,z) for a = 1 and complex b & z by inverse
  // evaluation of continued fractions
  const int Nmax = 1e3; // actually must be calculated from error estimation!

  complex_double b = CMPLX(*b_re - 1.0, *b_im), z = CMPLX(*z_re, *z_im);

  complex_double term;

  int n = Nmax;

  term = n * z / (b - z + n);

  for (n = Nmax - 1; n > 0; n--) {
    term = n * z / (b - z + n + term);
  }

  term = b / (b - z + term);

  *f_re = creal(term);
  *f_im = cimag(term);
}

/*******************************************************************/

void hypergeometric1f1_cont_fract_1_dir(double *b_re, double *b_im,
                                        double *z_re, double *z_im,
                                        double *f_re, double *f_im)
{
  // computes function 1F1(a,b,z) for a = 1 and complex b & z
  // by direct evaluation of continued fractions:
  // warning: this method is UNSTABLE sometimes!
  const int Nmax = 1e6;
  const double eps = DBL_EPSILON;
  double err;

  complex_double b = CMPLX(*b_re - 1.0, *b_im), z = CMPLX(*z_re, *z_im);

  complex_double An2, An1, An, Bn2, Bn1, Bn;
  complex_double an, bn;
  complex_double Sn1, Sn;

  An2 = 1.0e0;
  An1 = 0.0e0;
  Bn2 = 0.0e0;
  Bn1 = 1.0e0;
  Sn1 = An1 / Bn1;

  int n;

  for (n = 1; n < Nmax; n++) {
    an = n * z;
    bn = b - z + n;

    An = bn * An1 + an * An2;
    Bn = bn * Bn1 + an * Bn2;
    Sn = An / Bn;

    err = min(cabs((Sn - Sn1) / Sn), cabs(Sn - Sn1));
    if (err < eps) {
      break;
    }

    An2 = An1;
    An1 = An;
    Bn2 = Bn1;
    Bn1 = Bn;
    Sn1 = Sn;
  }

  Sn = b / (b - z + Sn);

  if (n == Nmax) {
    fprintf(stderr,
            "\nwarning: hypergeometric1f1_cont_fract_1_dir: n = %d err = "
            "%le\nSn = %le %le",
            n, err, creal(Sn), cimag(Sn));
    fprintf(stderr, "\nb = %le %le z = %le %le",
            creal(b) + 1.0, cimag(b), creal(z), cimag(z));
    fprintf(stderr, "\nBn = %le %le An = %le %le",
            creal(Bn), cimag(Bn), creal(An), cimag(An));
  }

  *f_re = creal(Sn);
  *f_im = cimag(Sn);
}

/*******************************************************************/
