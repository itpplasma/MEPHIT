#ifndef HYPER1F1_H
#define HYPER1F1_H

#ifdef __cplusplus
  #if __cplusplus >= 201103L
    #include <complex>
    typedef std::complex<double> complex_double;
  #else
    #error "need at least C++11 for binary compatibility with C99 _Complex"
  #endif
#else
  #if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #include <complex.h>
    typedef _Complex double complex_double;
  #else
    #error "need at least C99 for complex number support"
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

void hypergeometric1f1_quad(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_kummer_nmax(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_kummer_ada(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_kummer_modified_0_nmax(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_kummer_modified_0_ada(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_cont_fract_1_modified_0_ada(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_kummer_modified_0_accel(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_kummer_modified_1(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_cont_fract_1_dir(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_cont_fract_1_inv_nmax(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

void hypergeometric1f1_cont_fract_1_inv_ada(double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

#ifdef __cplusplus
}
#endif

#endif  // HYPER1F1_H
