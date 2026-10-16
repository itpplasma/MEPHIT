#include "mephit_native_mfem.h"
#include "miscellaneous.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
  mfem::SparseMatrix matrix(3, 3);
  matrix.Add(0, 0, 4.0); matrix.Add(0, 2, -1.0);
  matrix.Add(1, 0, 2.0); matrix.Add(1, 1, 3.0);
  matrix.Add(2, 1, -2.0); matrix.Add(2, 2, 5.0);
  matrix.Finalize();
  mfem::Vector rhs(3), solution(3);
  // Independent scalar equations for x=(1,-2,.5), not CSR conversion code.
  rhs(0) = 3.5; rhs(1) = -4.0; rhs(2) = 6.5;
  mephit::solve_native(matrix, rhs, solution);
  if (std::abs(solution(0)-1.0)>2e-14 || std::abs(solution(1)+2.0)>2e-14 ||
      std::abs(solution(2)-.5)>2e-14) { return 1; }
  mfem::Array<int> essential(1); essential[0]=1;
  mfem::Vector imposed(3); imposed=0.0; imposed(1)=-2.0;
  solution=utility::mfem_helper::solve_penalty(matrix,rhs,essential,imposed,8.0);
  if (std::abs(solution(0)-1.0)>2e-14 || std::abs(solution(1)+2.0)>2e-14 ||
      std::abs(solution(2)-.5)>2e-14) { return 4; }
  mfem::SparseMatrix singular(2, 2);
  singular.Add(0,0,1.0); singular.Add(1,0,2.0); singular.Finalize();
  mfem::Vector b(2), out(2); b=1.0; out=17.0;
  bool failed=false;
  try { mephit::solve_native(singular,b,out); }
  catch (const std::runtime_error &) { failed=true; }
  if (!failed || out(0)!=17.0 || out(1)!=17.0) { return 2; }
  failed=false;
  mfem::SparseMatrix unfinalized(2,2);
  try { mephit::solve_native(unfinalized,b,out); }
  catch (const std::runtime_error &) { failed=true; }
  if (!failed) { return 3; }
  std::cout << "PASS MFEM native CSR conversion, solve, singular and state checks\n";
}
