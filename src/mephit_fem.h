#ifndef MEPHIT_FEM_H
#define MEPHIT_FEM_H

#ifdef __cplusplus
  #include <complex>
  typedef complex<double> complex_double;
#else
  #if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #include <complex.h>
    typedef complex double complex_double;
  #else
    #error "need at least C99 for complex number support"
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

void FEM_init(const int tormode, const int nedge, const int runmode);
void FEM_extend_mesh(void);
void FEM_compute_Bn(const int nedge, const complex_double *restrict Jn, complex_double *restrict Bn);
void FEM_compute_L2int(const int nedge, const complex_double *restrict elem, double *restrict L2int);
void FEM_deinit(void);

#ifdef __cplusplus
}
#endif

#endif  // MEPHIT_FEM_H
