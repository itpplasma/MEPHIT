#ifndef MEPHIT_MATH_H
#define MEPHIT_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

// C bindings of the libneo math kit wrappers in math_c_api.f90.

void mephit_integrate_gk(double (*f)(double x, void *params), void *params,
                         double a, double b, double epsabs, double epsrel,
                         int limit, double *result, double *abserr, int *ierr);

void mephit_gauss_legendre(int n, double a, double b, double *x, double *w);

#ifdef __cplusplus
}
#endif

#endif  // MEPHIT_MATH_H
