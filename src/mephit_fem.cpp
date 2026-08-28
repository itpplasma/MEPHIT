#include "mephit_fem.h"
#include "mfem.hpp"
#ifdef USE_MFEM_MDE
#include "magnetic_differential_equation.h"
#endif  // USE_MFEM_MDE
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cstdio>
#include <vector>
#include <map>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::d2::point_xy<double> RZ_point;
typedef bg::model::box<RZ_point> RZ_box;
typedef std::pair<RZ_box, int> RZ_index;

static bgi::rtree< RZ_index, bgi::rstar<16> > RZ_tree;
static std::vector<int> query_results;

extern "C" void Rtree_init(int ntri, double *tri_bb)
{
  RZ_tree.clear();
  for (int k = 0; k < ntri; ++k) {
    // construct bounding box object
    RZ_box bb = RZ_box(RZ_point(tri_bb[4 * k], tri_bb[4 * k + 1]),
                       RZ_point(tri_bb[4 * k + 2], tri_bb[4 * k + 3]));
    RZ_tree.insert(std::make_pair(bb, k + 1));
  }
}

extern "C" void Rtree_query(double R, double Z, int *result_size, int **results)
{
  query_results.clear();
  for (auto it = RZ_tree.qbegin(bgi::contains(RZ_point(R, Z))); it != RZ_tree.qend(); ++it) {
    query_results.push_back(it->second);
  }
  *result_size = static_cast<int>(query_results.size());
  *results = query_results.data();
}

#ifdef USE_MFEM_MDE

typedef std::map<std::pair<double, double>, size_t> points_2D;

/*
void count_points_2D(points_2D& assoc, const double R, const double Z)
{
  auto const key = std::make_pair(R, Z);
  auto it = assoc.find(key);
  if (assoc.end() == it) {
    assoc.insert(it, std::make_pair(key, static_cast<size_t>(1)));
  } else {
    it->second += 1;
  }
}
*/

extern "C" int FEM_test(const char *mesh_file,
                        const int tor_mode,
                        const int n_dof,
                        complex_double *dof,
                        real_vector_field *unit_B0,
                        complex_scalar_field *MDE_inhom)
{
  int const fe_space_order = 1;
  int const quadrature_order = 2;
  try {
    points_2D points_lhs, points_rhs;
    // construct functors for real and imaginary part
    // (due to MFEM reasons, it has to be evaluated twice)
    auto const f_r = [MDE_inhom, &points_lhs](const double R, const double Z) {
      double f[2];
      MDE_inhom(R, Z, reinterpret_cast<complex_double *>(f));
      /* count_points_2D(points_lhs, R, Z); */
      return f[0];
    };
    auto const f_i = [MDE_inhom, &points_lhs](const double R, const double Z) {
      double f[2];
      MDE_inhom(R, Z, reinterpret_cast<complex_double *>(f));
      /* count_points_2D(points_lhs, R, Z); */
      return f[1];
    };
    auto const h_phi = [unit_B0, &points_rhs](const double R, const double Z) {
      double h[3];
      unit_B0(R, Z, h);
      /* count_points_2D(points_rhs, R, Z); */
      return h[2] / R;
    };
    auto const h_t = [unit_B0, &points_rhs](const double R, const double Z) {
      double h[3];
      unit_B0(R, Z, h);
      /* count_points_2D(points_rhs, R, Z); */
      mfem::Vector h_vec(2);
      h_vec(0) = h[0];
      h_vec(1) = h[1];
      return h_vec;
    };
    // read mesh
    mfem::Mesh mesh(mesh_file);
    // construct matrix and rhs
    // get system and rhs
    const int dim = mesh.Dimension();
    mfem::H1_FECollection u_coll(fe_space_order, dim);
    mfem::FiniteElementSpace U(&mesh, &u_coll);
    mfem::SparseMatrix Sys = mde::construct_mde_matrix(U, h_t, h_phi, tor_mode, quadrature_order);
    mfem::Vector rhs = mde::construct_mde_rhs(U, h_t, h_phi, f_r, f_i, tor_mode, quadrature_order);
    // solve system
    // construct solver
    // If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
    mfem::UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(Sys);
    mfem::Vector sol(rhs.Size());
    umf_solver.Mult(rhs, sol);
    if (n_dof != U.GetVSize()) {
      const int dof_mesh = U.GetVSize();
      fprintf(stderr, "Mismatch of vector sizes: n_dof = %i, GetVSize() -> %i.\n", n_dof, dof_mesh);
      return 1;
    }
    for (int i = 0; i < n_dof; ++i){
      reinterpret_cast<double *>(dof)[2 * i] = sol(i);
      reinterpret_cast<double *>(dof)[2 * i + 1] = sol(n_dof + i);
    }
    /*
    FILE *file;
    file = fopen("lhs.txt", "w");
    if (file != nullptr) {
      for (auto const& [key, value] : points_lhs) {
        fprintf(file, "%.16e %.16e %lu\n", key.first, key.second, value);
      }
      fclose(file);
    }
    file = fopen("rhs.txt", "w");
    if (file != nullptr) {
      for (auto const& [key, value] : points_rhs) {
        fprintf(file, "%.16e %.16e %lu\n", key.first, key.second, value);
      }
      fclose(file);
    }
    */
  } catch (...) {
    return 1;
  }
  return 0;
}

#endif  // USE_MFEM_MDE

/* made with help from ChatGPT */

class FourierGaugedCurlInterpolator : public mfem::DiscreteInterpolator
{
private:
  mfem::real_t n;

public:
  // factor should be -n for
  //
  // B_R = -n A_Z
  // B_Z =  n A_R
  //
  FourierGaugedCurlInterpolator(mfem::real_t tor_mode) : n(-tor_mode) {}

  void AssembleElementMatrix2(
    const mfem::FiniteElement &dom_fe,
    const mfem::FiniteElement &ran_fe,
    mfem::ElementTransformation &Trans,
    mfem::DenseMatrix &elmat) override
  {
    MFEM_VERIFY(dom_fe.GetDim() == 2 &&
                ran_fe.GetDim() == 2,
                "FourierGaugedCurlInterpolator: 2D elements required.");

    MFEM_VERIFY(dom_fe.GetMapType() == mfem::FiniteElement::H_CURL,
                "FourierGaugedCurlInterpolator: domain must be H(curl).");

    MFEM_VERIFY(ran_fe.GetMapType() == mfem::FiniteElement::H_DIV,
                "FourierGaugedCurlInterpolator: range must be H(div).");

    MFEM_VERIFY(dom_fe.GetDof() == ran_fe.GetDof(),
                "FourierGaugedCurlInterpolator: incompatible number of DOFs.");

    // Lowest-order ND_1 <-> RT_0:
    MFEM_VERIFY(dom_fe.GetDof() == 3,
                "FourierGaugedCurlInterpolator: "
                "this implementation assumes lowest-order elements.");

    elmat.SetSize(ran_fe.GetDof(), dom_fe.GetDof());
    elmat = 0.0;
    for (int i = 0; i < dom_fe.GetDof(); i++) {
      elmat(i, i) = n;
    }
  }
};

class MaxwellSolver {
public:
  const double c = 29979245800.0;
  // current formalism only works in lowest order
  const int order = 0;
  const int n;
  mfem::Mesh mesh;
  mfem::RT_FECollection RT0;
  mfem::ND_FECollection ND1;
  mfem::FiniteElementSpace Hdiv;
  mfem::FiniteElementSpace Hrot;
  mfem::BilinearForm potential;
  mfem::LinearForm source;
  mfem::DiscreteLinearOperator rot;
  mfem::Array<int> ess_tdof_list;
  mfem::GridFunction Hdiv_elem;
  mfem::VectorGridFunctionCoefficient Jn_interp;
  mfem::GridFunction An;
  mfem::Vector solution;
  mfem::Vector rhs;
  mfem::OperatorPtr lhs;
  mfem::UMFPackSolver umf;
  std::vector<int> edge_map;
  std::vector<int> sign_map;

  MaxwellSolver(const char* mesh_file, const int tor_mode);
  void map_edges(const char* edgemap_file);
  void assemble();
  void compute_magfn(const int nedge, const complex_double* Jn, complex_double* Bn);
};

MaxwellSolver::MaxwellSolver(const char* mesh_file, const int tor_mode)
  : n(tor_mode)
  // generate_edges = 0, refine = 0, fix_orientation = true
  // without refinement, local vertex (and thus edge) order is retained
  , mesh(mesh_file, 0, 0, true)
  , RT0(order, 2)
  , ND1(order + 1, 2)
  , Hdiv(&mesh, &RT0)
  , Hrot(&mesh, &ND1)
  , potential(&Hrot)
  , source(&Hrot)
  , rot(&Hrot, &Hdiv)
  , Hdiv_elem(&Hdiv)
  , Jn_interp(&Hdiv_elem)
  , An(&Hrot)
{}

void MaxwellSolver::map_edges(const char* edgemap_file)
{
  FILE* file;
  int ktri, ke, result, nedge;
  std::vector<int> mephit_ktri, mephit_ke;
  mfem::Array<int> edges, orient;
  file = fopen(edgemap_file, "r");
  nedge = 0;
  while (!feof(file)) {
    result = fscanf(file, "%d %d", &ktri, &ke);
    if (result != 2) break;
    nedge++;
    mephit_ktri.push_back(ktri - 1);
    mephit_ke.push_back(abs(ke) - 1);
  }
  fclose(file);
  edge_map.resize(nedge),
  sign_map.resize(nedge);
  for (int kedge = 0; kedge < nedge; kedge++) {
    mesh.GetElementEdges(mephit_ktri[kedge], edges, orient);
    edge_map[kedge] = edges[mephit_ke[kedge]];
    sign_map[kedge] = orient[mephit_ke[kedge]];
  }
}

void MaxwellSolver::assemble()
{
  mfem::FunctionCoefficient R(
    [](const mfem::Vector &X)
    {
      return X(0);
    }
  );
  mfem::CurlCurlIntegrator* const transverse_curl = new mfem::CurlCurlIntegrator(R);
  mfem::FunctionCoefficient n_squared_over_R(
    [this](const mfem::Vector &X)
    {
      return this->n * this->n / X(0);
    }
  );
  mfem::VectorFEMassIntegrator* const longitudinal_curl = new mfem::VectorFEMassIntegrator(n_squared_over_R);
  potential.AddDomainIntegrator(transverse_curl);
  potential.AddDomainIntegrator(longitudinal_curl);
  potential.Assemble();
  potential.Finalize();
  mfem::VectorFEDomainLFIntegrator* const curr_dens = new mfem::VectorFEDomainLFIntegrator(Jn_interp);
  source.AddDomainIntegrator(curr_dens);
  rot.AddDomainInterpolator(new FourierGaugedCurlInterpolator(n));
  rot.Assemble();
  rot.Finalize();
  mfem::Array<int> ess_bdr(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  Hrot.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  potential.FormSystemMatrix(ess_tdof_list, lhs);
}

void MaxwellSolver::compute_magfn(const int nedge, const complex_double* Jn, complex_double* Bn)
{
  for (ptrdiff_t im = 0; im <= 1; im++) {
    Hdiv_elem = 0.0;
    for (size_t k = 0; k < nedge; k++) {
      Hdiv_elem(edge_map[k]) = -0.25 * M_PI / c * sign_map[k] *
        reinterpret_cast<const double*>(Jn)[2 * k + im];
    }
    source.Assemble();
    An = 0.0;
    potential.FormLinearSystem(ess_tdof_list, An, source, lhs, solution, rhs);
    umf.SetOperator(dynamic_cast<mfem::SparseMatrix&>(*lhs));
    umf.Mult(rhs, solution);
    potential.RecoverFEMSolution(solution, source, An);
    rot.Mult(An, Hdiv_elem);
    for (size_t k = 0; k < nedge; k++) {
      // multiply by imaginary unit
      // im == 0:  Im B_n  <-  Re A_n
      // im == 1:  Re B_n  <- -Im A_n
      reinterpret_cast<double*>(Bn)[2 * k + (1 - im)] = (1 - 2 * im) *
        sign_map[k] * Hdiv_elem(edge_map[k]);
    }
  }
}

extern "C" void* MFEM_init(const int tor_mode, const char* mesh_file, const char* edgemap_file)
{
  MaxwellSolver* const maxwell_solver = new MaxwellSolver(mesh_file, tor_mode);
  maxwell_solver->map_edges(edgemap_file);
  maxwell_solver->assemble();
  return static_cast<void*>(maxwell_solver);
}

extern "C" void MFEM_compute_magfn(void* maxwell_solver, const int nedge, const complex_double* Jn, complex_double* Bn)
{
  if (maxwell_solver) {
    static_cast<MaxwellSolver*>(maxwell_solver)->compute_magfn(nedge, Jn, Bn);
  }
  return;
}

extern "C" void MFEM_deinit(void* maxwell_solver)
{
  if (maxwell_solver) {
    delete static_cast<MaxwellSolver*>(maxwell_solver);
  }
  return;
}
