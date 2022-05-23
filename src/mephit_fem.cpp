#include "mephit_fem.h"
#include <cstdio>

extern "C" int FEM_test(const char *mesh_file,
                        const int tor_mode,
                        const int n_dof,
                        complex_double *dof,
                        real_vector_field *unit_B0,
                        complex_scalar_field *MDE_inhom)
{
  double R, Z, h[3], f[2];
  try {
    printf("FEM_test:\n"
           "mesh_file = \"%s\"\n"
           "tor_mode = %i\n"
           "n_dof = %i\n",
           mesh_file, tor_mode, n_dof);
    R = 169.0; Z = 0.0;
    unit_B0(R, Z, h);
    MDE_inhom(R, Z, reinterpret_cast<complex_double *>(f));
    printf("R = %f, Z = %f, h = [%f, %f, %f], f = [%e, %e]\n",
           R, Z, h[0], h[1], h[2], f[0], f[1]);
    // (h_1 * partial_R + h_2 * partial_Z + 1i * tor_mode * h_3 / R) * u = f
    for (int k = 0; k < n_dof; ++k) {
      reinterpret_cast<double *>(dof)[2 * k] = static_cast<double>(k+1);
      reinterpret_cast<double *>(dof)[2 * k + 1] = static_cast<double>(-k-1);
    }
  } catch (...) {
    return 1;
  }
  return 0;
}
