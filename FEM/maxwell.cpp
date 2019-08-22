#include <mfem.hpp>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

const int dim = 2;
const int order = 1;
const int n = 2; // harmonic number

double f_r(const Vector & x)
{
   return x(0);
}

double f_nsq_over_r(const Vector & x)
{
   return n*1.0/x(0);
}

int main(int argc, char *argv[])
{
   bool doplot = true;
   cout << "Starting FEM" << endl;
   const char *mesh_file = "../PRELOAD/inputformaxwell_ext.mesh";
   const char *curr_file = "../MHD/currn.dat";

   Mesh *mesh = new Mesh(mesh_file, 1, 1);

   FiniteElementCollection *HcurlColl = new ND_FECollection(order, dim);
   FiniteElementSpace *Hcurl = new FiniteElementSpace(mesh, HcurlColl);

   GridFunction currRe(Hcurl);
   GridFunction currIm(Hcurl);

   currRe = 0.0;
   currIm = 0.0;

   ifstream currf;
   currf.open(curr_file, ios_base::in);
   int nt = mesh->GetNE();    // number of triangles
   for(int kt=0; kt<nt; kt++) 
   {
      double curr[8];
      Array<int> edges;
      Array<int> orient;
      mesh->GetElementEdges(kt, edges, orient);
      if (mesh->GetElement(kt)->GetAttribute() == 1)
      { 
         // read internal edge currents
         currf >> curr[0] >> curr[1] >> curr[2] >> curr[3]
               >> curr[4] >> curr[5] >> curr[6] >> curr[7];

         currRe[edges[0]] = curr[1]*orient[0];
         currRe[edges[1]] = curr[2]*orient[1];
         currRe[edges[2]] = curr[0]*orient[2];
      }
   }
   currf.close(); 

   Array<int> ess_tdof_list;
   if (mesh->bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 1;
      Hcurl->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }
/*
   // 6. Set up the linear form b(.) which corresponds to the right-hand side
   //    of the FEM linear system, which in this case is (f,phi_i) where f is
   //    given by the function f_exact and phi_i are the basis functions in the
   //    finite element fespace.
   LinearForm *b = new LinearForm(Hcurl);
   b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(currRe));
   b->Assemble();

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x by projecting the exact
   //    solution. Note that only values from the boundary edges will be used
   //    when eliminating the non-homogeneous boundary condition to modify the
   //    r.h.s. vector b.
   GridFunction x(fespace);
   VectorFunctionCoefficient E(sdim, E_exact);
   x.ProjectCoefficient(E);

   Coefficient *r = new FunctionCoefficient(f_r);
   Coefficient *nsq_over_r = new FunctionCoefficient(f_nsq_over_r);
   BilinearForm *a = new BilinearForm(Hcurl);
   a->AddDomainIntegrator(new CurlCurlIntegrator(*r));
   a->AddDomainIntegrator(new VectorFEMassIntegrator(*nsq_over_r));
   a->Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
*/
   if (doplot)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << currRe << flush;
   }
   
   delete mesh;

   return 0;
}
