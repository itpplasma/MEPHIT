#ifndef MEPHIT_FEM_H
#define MEPHIT_FEM_H

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

void FEM_init(const int tormode, const int nedge, const int npoint, const int runmode);
void FEM_extend_mesh(void);
void FEM_compute_magfn(const int nedge,
                       const int npoint,
                       const complex_double *Jn,
                       complex_double *Bn,
                       complex_double *AnR,
                       complex_double *AnZ);
void FEM_compute_L2int(const int nedge, const complex_double *elem, double *L2int);
void FEM_deinit(void);

void gauss_legendre_unit_interval(int order, double *points, double *weights);

void FEM_triangulate_external(const int npt_inner,
                              const int npt_outer,
                              const double *bdry_R,
                              const double *bdry_Z,
                              const double R_mid,
                              const double Z_mid,
                              const char *fname);
void Rtree_init(int ntri, double *tri_bb);
void Rtree_query(double R, double Z, int *result_size, int **results);

typedef void real_vector_field(const double R,
                               const double Z,
                               double *vector);
typedef void complex_scalar_field(const double R,
                                  const double Z,
                                  complex_double *scalar);

#ifdef USE_MFEM
int FEM_test(const char *mesh_file,
             const int tor_mode,
             const int n_dof,
             complex_double *dof,
             real_vector_field *unit_B0,
             complex_scalar_field *MDE_inhom);
#endif  // USE_MFEM

#ifdef __cplusplus
}
#endif

#endif  // MEPHIT_FEM_H
