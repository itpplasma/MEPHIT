//
// Created by michi on 02.02.22.
//

#include "miscellaneous.h"
#include <random>

double utility::mfem_helper::evaluate_bilinear_form(const mfem::SparseMatrix &A, const mfem::Vector &u,
                                                    const mfem::Vector &v) {

    const int n_dom = A.Width();
    const int n_cod = A.Height();
    if(n_dom != u.Size()){
        throw std::runtime_error("Size of u does not match domain");
    }
    if(n_cod != v.Size()){
        throw std::runtime_error("Size of v does not match codomain");
    }
    mfem::Vector t(n_cod);
    A.Mult(u, t);
    const double res = t.operator*(v);
    return res;
}

mfem::Vector utility::mfem_helper::get_random_vector(const int n) {
    mfem::Vector rv(n);
    //randomly populate vector
    double lower_bound = 0.;
    double upper_bound = 10.;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;
    for(int i = 0; i<n; ++i){
        rv(i) = unif(re);
    }
    return rv;
}

double utility::mfem_helper::get_operator_application_difference(mfem::Operator *A, mfem::Operator *B) {
    if(A->Height() != B->Height()){
        int const height_a = A->Height();
        int const height_b = B->Height();
        throw std::runtime_error("get_operator_application_difference: Operators must have same height!");
    }
    if(A->Width() != B->Width()){
        int const width_a = A->Width();
        int const width_b = B->Width();
        throw std::runtime_error("get_operator_application_difference: Operators must have same width!");
    }
    const int n = A->Width();
    const int m = A->Height();
    mfem::Vector rv = get_random_vector(n);
    mfem::Vector res1(m), res2(m);
    A->Mult(rv, res1);
    B->Mult(rv, res2);

    res1 -= res2;

    const double l2_error = res1.Norml2();

    return l2_error;
}

double utility::mfem_helper::evaluate_bilinear_form(mfem::Operator *A, const mfem::Vector &u, const mfem::Vector &v) {

    const int n_dom = A->Width();
    const int n_cod = A->Height();
    if(n_dom != u.Size()){
        throw std::runtime_error("Size of u does not match domain");
    }
    if(n_cod != v.Size()){
        throw std::runtime_error("Size of v does not match codomain");
    }
    mfem::Vector t(n_cod);
    A->Mult(u, t);
    const double res = t.operator*(v);
    return res;
}

mfem::SparseMatrix utility::mfem_helper::block_to_sparse_matrix(const mfem::BlockMatrix &B) {
    //more or less directly copied from BlockMatrix::PrintMatlab
    mfem::SparseMatrix S(B.Height(), B.Width());
    mfem::Vector row_data;
    mfem::Array<int> row_ind;
    int nnz_elem = B.NumNonZeroElems();
    int i, j;
    for (i = 0; i < B.RowOffsets().Last(); i++)
    {
        B.GetRow(i, row_ind, row_data);
        for (j = 0; j < row_ind.Size(); j++)
        {
            S.Add(i, row_ind[j], row_data[j]);
        }
    }
    //S.Finalize();
    return S;
}

mfem::Vector utility::mfem_helper::solve_penalty(const mfem::SparseMatrix &A, const mfem::Vector &rhs,
                                                 const mfem::Array<int> &essential_dofs, const mfem::Vector &u,
                                                 const double penalty_fac) {

    mfem::SparseMatrix P(A);
    mfem::Vector b(rhs);
    for(int i=0; i<essential_dofs.Size(); ++i){
        const int dof = essential_dofs[i];
        P.Add(dof,dof, penalty_fac);
        b(dof) += u(dof)*penalty_fac;
    }
    P.Finalize();
    mfem::Vector sol(rhs.Size());

    //construct solver
#ifndef MFEM_USE_SUITESPARSE
    throw std::runtime_error("solver not available. Recopile mfem!");
#else
    // If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
    mfem::UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(P);
    umf_solver.Mult(b, sol);
#endif

    return sol;
}

